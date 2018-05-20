# vim: sw=2 ts=2 sts=2 tw=80 et:
# Nim re-write of FALCON_unzip/falcon_unzip/rr_hctg_track.py
import "../msgpack4nim/msgpack4nim.nim"
import "../msgpack4nim/msgpack4collection.nim"
from "../nim-heap/binaryheap" import nil
from sys import log

from algorithm import nil
from hashes import nil
from json import nil
import json
from os import nil
from ospaths import nil
from osproc import nil
#from "pure/osproc.nim" import nil
from parsecsv import nil
from strutils import splitWhitespace, strip, `%`
from sequtils import nil
from sets import incl, contains, items
from streams import write, setPosition, readAll, close, writeLine
import sequtils
import tables
import threadpool

#import os, tables, strutils, streams, algorithm, json, hashes
#import memfiles


type
  # A ctg is a 9-char array.
  str9 = array[0..9, char]
  #mytuple = tuple[score: int, rid: string]
  mytuple = tuple[score: int, rid: str9] #overlap_len, q_id
  mytable = tables.Table[string, seq[mytuple]]
  myprioritytable = tables.Table[string, binaryheap.Heap[mytuple]]

proc toseq(src: str9): seq[int] =
  var dst = newSeq[int](9)
  assert len(dst) == 9
  for i in 0..<9:
    dst[i] = src[i].int
  return dst
proc fromstring(src: string): str9 =
  # dst must be same length as src
  assert len(src) == 9
  for i in 0..<9:
    result[i] = src[i]
  #log("fromstring('", src, "'->'", $toseq(result))

proc tostring(src: str9): string =
  var dst = newString(9)
  for i in 0..<9:
    dst[i] = src[i]
  return dst

proc tostring(src: str9, dst: var string) =
  assert len(dst) == 9
  for i in 0..<9:
    dst[i] = src[i]

proc hash(x: str9): hashes.Hash =
  return hashes.hash(tostring(x))

proc unpack_type*(ss: var msgpack4nim.MsgStream, x: var str9) =
  var str: string
  msgpack4nim.unpack(ss, str)
  copyMem(addr x[0], addr str[0], 9)

proc pack_type*(ss: var msgpack4nim.MsgStream, x: str9) =
  let str: string = tostring(x)
  msgpack4nim.pack(ss, str)

proc pack_type*(ss: var msgpack4nim.MsgStream, x: binaryheap.Heap[mytuple]) =
  let xseq: seq[mytuple] = sequtils.toSeq(binaryheap.items(x)) # unsorted
  msgpack4nim.pack(ss, xseq)
proc `$`(x: str9): string =
  return tostring(x)
proc `$`(x: mytuple): string =
  return "[$#, $#]" % [$x.score, $x.rid]
#[
proc `$`(x: Table): string =
  var ss = streams.newStringStream()
  ss.write("len=" & $len(x))
  ss.write("{")
  var first = true
  for k,v in items(x):
    if first:
      first = false
    else:
      ss.write(",")
    ss.write("\n\t'")
    ss.write($k)
    ss.write("': ")
    ss.write("<")
    ss.write($v)
    ss.write(">")
  ss.write("\n}")
  ss.setPosition(0)
  return ss.readAll()
]#
proc `$`(x: binaryheap.Heap[mytuple]): string =
  var ss = streams.newStringStream()
  ss.write("len=" & $binaryheap.size(x))
  ss.write("[")
  var first = true
  for v in binaryheap.items(x):
    if first:
      first = false
    else:
      ss.write(", ")
    ss.write("'")
    ss.write($v)
    ss.write("'")
  ss.write("]")
  ss.setPosition(0)
  return ss.readAll()


proc mytuple_cmp(x, y: mytuple): int =
  if x.score < y.score: return -1
  if x.score > y.score: return 1
  for i in 0..<9:
    if x.rid[i] != y.rid[i]: return system.cmp(x.rid[i], y.rid[i])
  return 0


proc glob_raw_reads_las*(): seq[string] =
  result = @[]
  # TODO: Use 0-rawreads/las-gather/las_fofn.json, which are all called 'merged.las' now.
  for fn in os.walkFiles("0-rawreads/las-merge-runs/m_*/uow-*/merged.las"):
    result.add(fn)
  algorithm.sort(result, system.cmp) # probably not needed

proc get_fns_from_json*(json_fn: string): seq[string] =
  log("Reading JSON from '", json_fn, "'")
  let jcontent = system.readFile(json_fn)
  result = to(json.parseJson(jcontent), seq[string]) # json.to()
  let fofn_dir = ospaths.parentDir(json_fn) # "" is ok
  for i in 0 ..< result.len:
    let fn = result[i]
    if not ospaths.isAbsolute(fn):
      result[i] = ospaths.joinPath(fofn_dir, fn)

proc fillTuple[T: tuple, V](input: openarray[V]): T =
  #assert input.len == len(result.fields)
  var index = 0
  for field in result.fields:
    assert input.len > index
    field = input[index]
    inc(index)
  assert len(input) == index
proc stream_get_rid_to_ctg_slow(infile: streams.Stream): auto =
  log("In stream_get_rid_to_ctg_slow()")
  type TRow = tuple[pid, rid, oid, ctg: string]
  let
    #rid_to_ctg = tables.newTable[string, seq[string]]()
    rid_to_ctg = newTable[string, sets.HashSet[string]](1 shl 20)
  var line: TaintedString = ""
  while streams.readLine(infile, line):
    let s = line.splitWhitespace()
    if len(s) == 0: continue
    let t: TRow = fillTuple[TRow, string](s)
    #if not rid_to_ctg.hasKey(t.rid):
    #  rid_to_ctg[t.rid] = newSeq[string]()
    #rid_to_ctg[t.rid].add(t.ctg)
    if not rid_to_ctg.contains(t.rid):
      rid_to_ctg[t.rid] = sets.initSet[string](2)
    rid_to_ctg[t.rid].incl(t.ctg)
  return rid_to_ctg

var logged_rid_to_ctg = false

proc stream_get_rid_to_ctg(ss: streams.Stream, fn="rid_to_ctg"): auto =
  log("In stream_get_rid_to_ctg()")
  var
    csv: parsecsv.CsvParser
    rid_to_ctg = newTable[string, sets.HashSet[string]](1 shl 20)

  parsecsv.open(csv, ss, fn, separator=' ')

  while parsecsv.readRow(csv):
    #if csv.row.len() < 4:
    #  continue
    if csv.row.len() != 4:
      raise newException(ValueError, $csv.row)
    #let pid = csv.row[0]
    let rid: string = csv.row[1]
    #let oid = csv.row[2]
    let ctg = csv.row[3]
    #log("rid=", rid, ", ctg=", ctg)
    if not rid_to_ctg.contains(rid):
      rid_to_ctg[rid] = sets.initSet[string](2)
    if (not logged_rid_to_ctg) and (ctg in rid_to_ctg[rid]):
      logged_rid_to_ctg = true
      log("In ", $fn, ", found dup ctg ", $ctg, " for rid ", $rid, " from Row ", $csv.row, " matching one of ", $sequtils.toSeq(rid_to_ctg[rid].items()))
    #assert ctg notin rid_to_ctg[rid]
    rid_to_ctg[rid].incl(ctg)
  return rid_to_ctg

proc get_rid_to_ctg_slow(fn: string): auto =
  log("In get_rid_to_ctg_slow(): len=", sys.filelength(fn), " for file '", fn, "'")
  let
    ss = streams.newFileStream(fn)
  assert(not ss.isNil())
  defer:
    streams.close(ss)
  return stream_get_rid_to_ctg_slow(ss)

proc get_rid_to_ctg*(fn: string): auto =
  #let content = readFile(fn)
  log("In get_rid_to_ctg(): len=", sys.filelength(fn), " for file '", fn, "'")
  let
    ss = streams.newFileStream(fn, system.fmRead)
  defer:
    streams.close(ss)
  return stream_get_rid_to_ctg(ss, fn)

#[
  #let mfile = memfiles.open(fn)
  #for line in memfiles.lines(mfile):
  #for line in lines(fn):
  #  echo line
]#


#min_len, bestn, rid_to_ctg, rid_to_phase
type
  Settings* = object {.requiresInit.}
    n_core*: int #48
    file_list*: seq[string]
    db_fn*: string #"?/raw_reads.db"
    phased_read_file*: string # ./3-unzip/all_phased_reads
    read_to_contig_map*: string # ./4-quiver/read_maps/read_to_contig_map
    rawread_ids*: string # ./2-asm-falcon/read_maps/dump_rawread_ids/rawread_ids
    output*: string # ./2-asm-falcon/read_maps/dump_rawread_ids/rawread_to_contigs
    min_len*: int # 2500
    stream*: bool # false
    debug*: bool # false (also single-threaded)
    bestn*: int # 40
  Phase = tuple[
    ctg_id: string,
    blockn: int16,
    phase: int16,
  ] # sizeof == 8 + 2 + 2 + pad == 16

proc stream_get_oid_to_phase(ss: streams.Stream, fn: string): auto =
  let oid_to_phase = tables.newTable[string, Phase]()
  log("In stream_get_oid_to_phase()")
  var
    csv: parsecsv.CsvParser

  parsecsv.open(csv, ss, fn, separator=' ')

  while parsecsv.readRow(csv):
    if csv.row.len() == 0:
      continue
    let p: Phase = (
      ctg_id: csv.row[1],
      blockn: strutils.parseInt(csv.row[2]).int16,
      phase: strutils.parseInt(csv.row[3]).int16
    )
    let oid = csv.row[6]
    oid_to_phase[oid] = p
  return oid_to_phase

proc stream_get_oid_to_phase_slow(ss: streams.Stream, fn: string): auto =
  let oid_to_phase = tables.newTable[string, Phase]()
  log("In stream_get_oid_to_phase_slow()")

  var line: TaintedString = ""
  while streams.readLine(ss, line):
    let row = line.splitWhitespace()
    if len(row) == 0: continue
    #let t: TRow = fillTuple[TRow, string](row)

    #log("phased_reads_row:", $row)
    let p: Phase = (
      ctg_id: row[1],
      blockn: strutils.parseInt(row[2]).int16,
      phase: strutils.parseInt(row[3]).int16
    )
    let oid = row[6]
    oid_to_phase[oid] = p
  return oid_to_phase

proc stream_get_rid_to_phase(rawread_ids_fn: string,
    oid_to_phase: TableRef[string, Phase]): auto =
  log("In stream_get_rid_to_phase_slow()")
  #let rid_to_oid: seq[string] = sequtils.toSeq(open(rawread_ids_fn).readAll().splitLines())
  let rid_to_oid: seq[string] = sequtils.toSeq(lines(open(rawread_ids_fn))) # 5x faster
  log(" Read ", len(rid_to_oid), " strings into rid_to_oid list.")
  let n_rids = rid_to_oid.len
  var rid_to_phase: ref seq[Phase]
  new(rid_to_phase)
  newSeq[Phase](rid_to_phase[], n_rids)

  var nin = 0
  var nout = 0
  for rid in 0..<n_rids:
    let oid = rid_to_oid[rid]
    rid_to_phase[rid] = tables.getOrDefault(oid_to_phase, oid)
    # We could save memory in the table by using "ref Phase", but
    # then the Phases would be reference-counted. Not a good trade-off.
    # And the real memory problem comes from the threads/sub-procs.
    if rid_to_phase[rid].ctg_id.is_nil:
      nout.inc
    else:
      nin.inc
    #log(" r2p[", rid, "]=>", rid_to_phase[rid], " for oid=", oid)
  log(" Of ", n_rids, ", in=", nin, " & out=", nout)
  return rid_to_phase

proc get_rid_to_phase(phased_reads_fn, rawread_ids_fn: string): auto =
  let fn0 = phased_reads_fn
  log("Reading phased_reads file: len=", sys.filelength(fn0), " for file '", fn0, "'")
  let
    ss0 = streams.newFileStream(fn0, system.fmRead)
  defer:
    streams.close(ss0)
  log(" sizeof(Phase):", sizeof(Phase))
  let oid_to_phase = stream_get_oid_to_phase(ss0, fn0)
  log(" len(oid_to_phase)=", len(oid_to_phase))

  let fn1 = rawread_ids_fn
  log("Reading rawread_ids file: len=", sys.filelength(fn1), " for file '", fn1, "'")
  #let
  #  ss1 = streams.newFileStream(fn1, system.fmRead)
  #defer:
  #  streams.close(ss1)
  let rid_to_phase = stream_get_rid_to_phase(fn1, oid_to_phase)
  log(" len(rid_to_phase)=", len(rid_to_phase[]))
  return rid_to_phase

var global_rid_to_ctg: ref Table[string, sets.HashSet[string]]
var global_rid_to_phase: ref seq[Phase]

proc tr_stage1(la4falcon_stream: streams.Stream, fn: string, min_len, bestn: int): myprioritytable =
  # for each read in the b-read column inside the LAS files, we
  #keep top `bestn` hits with a priority queue through all overlaps
  result = tables.initTable[string, binaryheap.Heap[mytuple]](1024)

  var
    csv: parsecsv.CsvParser

  parsecsv.open(csv, la4falcon_stream, fn, separator=' ')

  while parsecsv.readRow(csv):
    #log("row:", $csv.row)
    let
        q_id = csv.row[0]
        t_id = csv.row[1]
        overlap_len = -strutils.parseInt(csv.row[2])
        idt = strutils.parseFloat(csv.row[3])
        (q_s, q_e, q_l) = (strutils.parseInt(csv.row[5]), strutils.parseInt(csv.row[6]), strutils.parseInt(csv.row[7]))
        (t_s, t_e, t_l) = (strutils.parseInt(csv.row[9]), strutils.parseInt(csv.row[10]), strutils.parseInt(csv.row[11]))
    if t_l < min_len:
        continue
    if q_id notin global_rid_to_ctg:
        continue
    let t_id_int = strutils.parseInt(t_id)
    #log(" t_id_int:", t_id_int, ", len(r2p):", len(global_rid_to_phase))
    if t_id_int < len(global_rid_to_phase[]):
      let t_phase = global_rid_to_phase[t_id_int]
      if t_phase.ctg_id != nil:
        if t_phase.blockn != -1:
            let q_id_int = strutils.parseInt(q_id)
            #log("  q_id_int:", q_id_int)
            if q_id_int < len(global_rid_to_phase[]):
              let q_phase = global_rid_to_phase[q_id_int]
              if q_phase.ctg_id != nil:
                #log("   t_phase:", t_phase)
                #log("   q_phase:", q_phase)
                if (q_phase.ctg_id == t_phase.ctg_id and q_phase.blockn == t_phase.blockn and
                    q_phase.phase != t_phase.phase):
                  continue
    if t_id notin result:
      result[t_id] = binaryheap.newHeap[mytuple](mytuple_cmp)
    let it: mytuple = (score: overlap_len, rid: fromstring(q_id))
    #log("  rid:", $fromstring(q_id))
    #log("  it:", $it)
    if binaryheap.size(result[t_id]) < bestn:
      binaryheap.push(result[t_id], it)
      #log("  heap:", $result[t_id])
    else:
      discard binaryheap.pushPop(result[t_id], it)
      #log("  heap now:", $result[t_id])

proc captureCmd(cmd: string): auto =
  #log("$$($#)" % [cmd])
  let options: set[osproc.ProcessOption] = {osproc.poUsePath, osproc.poEvalCommand, osproc.poEchoCmd}
  result = osproc.execProcess(cmd, options=options)
  log(" len(stdout)==", len(result))

proc run_tr_stage1(db_fn, fn: string, min_len, bestn: int): string =
  let cmd = "LA4Falcon -m $# $#" % [db_fn, fn]
  #let la4falcon_output = captureCmd(cmd)
  #log("Parsing LA4Falcon output of size=", len(la4falcon_output), " bytes")
  #let infile = streams.newStringStream(la4falcon_output)
  #defer: infile.close()
  log("Parsing LA4Falcon output, streamed ...")
  let options: set[osproc.ProcessOption] = {osproc.poUsePath, osproc.poEvalCommand, osproc.poEchoCmd}
  let la4falcon_proc = osproc.startProcess(cmd, options=options)
  let infile = osproc.outputStream(la4falcon_proc)
  let rtn = tr_stage1(infile, fn, min_len, bestn)
  let rc = osproc.peekExitCode(la4falcon_proc)
  if rc != -1:
    log("LA4Falcon exited with:", rc)
  else:
    log("LA4Falcon has not exited. Terminating.")
  osproc.close(la4falcon_proc)
  let fn_rtn = "$#.rr_hctg_track.partial.msgpack" % [fn]
  log("Serialize '$#'" % [fn_rtn])
  var msgss = msgpack4nim.initMsgStream() #you can pass some buffer capacity here https://github.com/jangko/msgpack4nim/issues/18
  msgpack4nim.pack(msgss, rtn)
  var ss = streams.newFileStream(fn_rtn, system.fmWrite)
  defer: ss.close()
  ss.write(msgss.data)
  return fn_rtn

proc spawned_tr_stage1(db_fn, fn: string, min_len, bestn: int): string {.thread} =
  log("In new spawned_tr_stage1:", fn)
  #os.sleep(1000)
  #log("Returning...")
  #return "hello"
  return run_tr_stage1(db_fn, fn, min_len, bestn)

proc run_track_reads*(settings: Settings) =
  let rid_to_ctg = get_rid_to_ctg(settings.read_to_contig_map)
  let rid_to_phase = get_rid_to_phase(settings.phased_read_file, settings.rawread_ids)
  global_rid_to_ctg = rid_to_ctg
  global_rid_to_phase = rid_to_phase
  #log("spinning...")
  #while true:
  #  discard
  #if true:
  #  log("Early return")
  #  return
  var nsubs = len(settings.file_list)
  var responsesB = newSeq[FlowVarBase](nsubs)
  var responsesT = newSeq[FlowVar[string]](nsubs)
  for i in 0..<nsubs:
    let fn = settings.file_list[i]
    let fv: FlowVar[string] = spawn spawned_tr_stage1(settings.db_fn, fn, settings.min_len, settings.bestn)
    responsesT[i] = fv
    responsesB[i] = fv
  sync()
  #while nsubs > 0:
  #  log("Awaiting...")
  #  let i = awaitAny(responsesB)
  #  log(" Got i=", i)
  #  #assert i != -1
  #  let fv: FlowVar[string] = responsesT[i]
  #  let fn_rtn = ^fv
  #  log(" Wrote ", sys.filelength(fn_rtn), " bytes into msgpack file:", fn_rtn)
  #  nsubs.dec
  for i in 0..<nsubs:
    let fv: FlowVar[string] = responsesT[i]
    let fn_rtn = ^fv
    log(i, ": Wrote ", sys.filelength(fn_rtn), " bytes into msgpack file:", fn_rtn)

proc run_stage2*(
    test=false,
    read_to_contig_map="./4-quiver/track_reads/read_to_contig_map",
    #las_fofn_fn="./0-rawreads/las-gather/las_fofn.json",
    partials_fn="./4-quiver/track-reads/partials.json",
    output="./4-quiver/track-reads/rawread_to_contigs",
    bestn=40,
    # All the rest are ignored....
    n_core=0, phased_read_file="", rawread_ids="", min_len=2500,
    stream=false, debug=false, silent=false) =
  ## Stage 2 of rr_hctg_track
  if test:
    log("no tests")
    system.quit(system.QuitSuccess)
  #var fn_rtns: seq[string] = @[]
  let fn_rtns: seq[string] = get_fns_from_json(partials_fn)
  #for fn in file_list:
  #  fn_rtns.add(fn & ".rr_hctg_track.partial.msgpack")
  #log("file_list:", $file_list)
  log("inputs:", $fn_rtns)

  #Aggregate hits from each individual LAS and keep the best n hit.
  #Note that this does not guarantee that the final results is globally the best n hits espcially
  #when the number of `bestn` is too small.  In those case, if there is more hits from single LAS
  #file, then we will miss some good  hits.
  var bread_to_areads: myprioritytable = tables.initTable[string, binaryheap.Heap[mytuple]](1024)
  for fn_rtn in fn_rtns:
    when false:
      log("Deserialize ", fn_rtn, " as json...")
      #var jp: json.JsonParser
      #let infile = streams.newFileStream(fn_rtn, system.fmRead)
      #defer:
      #  infile.close()
      #json.open(jp, infile, fn_rtn)
      let jcontent = system.readFile(fn_rtn)
      let j = json.parseJson(jcontent)
    when true:
      log("Deserialize ", fn_rtn, " as msgpack...")
      doAssert(os.fileExists(fn_rtn), fn_rtn)
      let infile = streams.newFileStream(fn_rtn, system.fmRead)
      defer:
        infile.close()
      var msgss = msgpack4nim.initMsgStream()
      msgss.data = infile.readAll()
      var res: mytable
      msgpack4nim.unpack(msgss, res)
      log("  Updating bread_to_areads table with ", len(res), " breads.")
      for k,v in res.pairs:
        if not bread_to_areads.contains(k):
          #echo k, " not in b2a"
          bread_to_areads[k] = binaryheap.newHeap[mytuple](mytuple_cmp)
          #echo "len(b2a)=", len(bread_to_areads)
        for item in v:
          if binaryheap.size(bread_to_areads[k]) < bestn:
            binaryheap.push(bread_to_areads[k], item)
          else:
            discard binaryheap.pushPop(bread_to_areads[k], item)
      log("  Num breads now:", len(bread_to_areads))

  # rid_to_oid can be helpful for debugging, but otherwise we do not need it.
  #rid_to_oid = open(os.path.join(rawread_dir, 'dump_rawread_ids', 'rawread_ids')).read().split('\n')

  let
    rid_to_ctg = get_rid_to_ctg(read_to_contig_map)
    rawread_to_contigs_fn = output
  log("Opening ", sys.abspath(rawread_to_contigs_fn), " for output.")
  let
    out_f = streams.newFileStream(rawread_to_contigs_fn, system.fmWrite)
  defer:
    out_f.close()

  type ScoreCount = tuple[score: int, count: int]
  proc `$`(val: ScoreCount): auto =
    return "(" & $val[0] & ", " & $val[1] & ")"

  log(" Writing ctg/scores for ", len(bread_to_areads), " breads.")
  var rid = newString(9) # same length as str9

  # For each b-read, we find the best contig map throgh the b->a->contig map.
  for bread in tables.keys(bread_to_areads):
    # Note that breads will be in arbitrary order, so the output will match
    # Python but in a different order.
    var ctg_score = tables.initOrderedTable[string, ScoreCount](1024)
    for x in binaryheap.items(bread_to_areads[bread]):
      tostring(x.rid, rid)
      if not rid_to_ctg.contains(rid):
        continue
      for ctg in rid_to_ctg[rid]:
        #discard tables.mgetOrPut(ctg_score, ctg, [0,0])
        if not ctg_score.contains(ctg):
          ctg_score[ctg] = (0, 0)
        ctg_score[ctg].score += -x.score
        ctg_score[ctg].count += 1
    ##oid = rid_to_oid[int(bread)]
    #type ctg_score_pair = tuple[ctg: string, score_count: array[2, int]]
    #var ctg_scores: seq[ctg_score_pair] = sequtils.toSeq(tables.pairs(ctg_score))
    #algorithm.sort(ctg_scores) do (x, y: ctg_score_pair) -> int:
    #  result = x.score_count[0]
    ctg_score.sort(proc (left, right: (string, ScoreCount)): int =
      if left[1].score != right[1].score:
        return system.cmp(left[1].score, right[1].score)
      else:
        return system.cmp(left[0], right[0])
    )
    var rank = 0
    for ctg, score_count in ctg_score.pairs:
      var in_ctg: int
      if rid_to_ctg.contains(bread) and rid_to_ctg[bread].contains(ctg):
        in_ctg = 1
      else:
        in_ctg = 0
      #out_f.writeLine(bread, " ", oid, " ", ctg, " ", score_count.count, " ", rank, " ", score_count.score, " ", in_ctg)
      out_f.writeLine(bread, " ", ctg, " ", score_count.count, " ", rank, " ", score_count.score, " ", in_ctg)
      rank += 1
  #log("spinning...")
  #while true:
  #  discard


proc test() =
  # First, test an assumption.
  assert 18602 == strutils.parseInt("000018602")

  # Test rid_to_ctg.
  let data_get_rid_to_ctg = """
000000026 000000026 ref1/27/0_39949 000000F
000000026 000000026 ref1/27/0_39949 000000C
000000278 000000287 ref2/128/0_42163 000000F_02
"""
  block:
    # slow
    let res = stream_get_rid_to_ctg_slow(streams.newStringStream(
        data_get_rid_to_ctg))
    assert len(res) == 2, $len(res)
    assert $res["000000026"] == "@[000000F, 000000C]", $res["000000026"]
    assert $res["000000287"] == "@[000000F_02]", $res["000000287"]
    echo "Passed."
  block:
    # faster, I think
    let res = stream_get_rid_to_ctg(streams.newStringStream(
        data_get_rid_to_ctg))
    assert len(res) == 2, $len(res)
    assert $res["000000026"] == "@[000000F, 000000C]", $res["000000026"]
    assert $res["000000287"] == "@[000000F_02]", $res["000000287"]
  echo "Success."
  block:
    let fn = "/lustre/hpcprod/cdunn/jira/se-809/run/4-quiver/read_maps/read_to_contig_map"
    #let res = get_rid_to_ctg(fn)

proc main() =
  test()

if isMainModule:
  test()
