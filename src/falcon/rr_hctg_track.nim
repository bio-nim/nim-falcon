# vim: sw=2 ts=2 sts=2 tw=80 et:
# Nim re-write of FALCON_unzip/falcon_unzip/rr_hctg_track.py
from "nim-heap/binaryheap" import nil
from algorithm import sort
from strutils import splitWhitespace, strip, `%`
from sequtils import nil
from streams import nil
import sequtils
import tables

proc fillTuple[T: tuple, V](input: openarray[V]): T =
  #assert input.len == len(result.fields)
  var index = 0
  for field in result.fields:
    assert input.len > index
    field = input[index]
    inc(index)
  assert len(input) == index
proc stream_get_rid_to_ctg(ss: streams.Stream): auto =
  type TRow = tuple[pid, rid, oid, ctg: string]
  let rid_to_ctg = tables.newTable[string, seq[string]]()
  var line: TaintedString = ""
  while streams.readLine(ss, line):
    let s = line.splitWhitespace()
    if len(s) == 0: continue
    let t: TRow = fillTuple[TRow, string](s)
    if not rid_to_ctg.hasKey(t.rid):
      rid_to_ctg[t.rid] = newSeq[string]()
    rid_to_ctg[t.rid].add(t.ctg)
  return rid_to_ctg
proc get_rid_to_ctg(fn: string): auto =
  let ss = streams.newFileStream(fn)
  assert(not ss.isNil())
  defer: streams.close(ss)
  return stream_get_rid_to_ctg(ss)
  
#min_len, bestn, rid_to_ctg, rid_to_phase
type
  Settings = object
    n_core: int #48
    fofn: string
    db_fn: string #"?/raw_reads.db"
    phased_read_file: string # ./3-unzip/all_phased_reads
    read_to_contig_map: string # ./4-quiver/read_maps/read_to_contig_map
    rawread_ids: string # ./2-asm-falcon/read_maps/dump_rawread_ids/rawread_ids
    output: string # ./2-asm-falcon/read_maps/dump_rawread_ids/rawread_to_contigs
    min_len: int # 2500
    stream: bool # false
    debug: bool # false (also single-threaded)
    bestn: int # 40
  Phase = object
    ctg_id: string
    blockn: int
    phase: int
  HeapItem = tuple[overlap_len: int, q_id: string] # overlap_len, q_id
  MyHeap = binaryheap.Heap[HeapItem]
  String2Heap = tables.Table[string, ref MyHeap]
proc cmp_heap_items(a, b: HeapItem): int =
  return 0 # a-b
proc run_track_reads(settings: Settings) =
  let rid_to_ctg = get_rid_to_ctg(settings.read_to_contig_map)
  let oid_to_phase = tables.newTable[string, ref Phase]()
  let f = open(settings.phased_read_file)
  assert(not f.isNil)
  defer: f.close()
  for line in lines(f):
    let row = line.strip().splitWhitespace()
    echo $row
    var p: ref Phase = new(Phase)
    p.ctg_id = row[1]
    p.blockn = strutils.parseInt(row[2])
    p.phase = strutils.parseInt(row[3])
    let oid = row[6]
    oid_to_phase[oid] = p

  let rid_to_oid = sequtils.toSeq(lines(open(settings.rawread_ids)))
  let n_rids = rid_to_oid.len
  var rid_to_phase = newSeq[ref Phase](n_rids)

  for rid in 0..<n_rids:
    let oid = rid_to_oid[rid]
    rid_to_phase[rid] = tables.getOrDefault(oid_to_phase, oid)
  #[
  Aggregate hits from each individual LAS and keep the best n hit.
  Note that this does not guarantee that the final results is globally the best n hits espcially
  when the number of `bestn` is too small.  In those case, if there is more hits from single LAS
  file, then we will miss some good  hits.
  ]#
  var res: ref String2Heap = tables.newTable[string, ref MyHeap]()
  let bread_to_areads: ref String2Heap = tables.newTable[string, ref MyHeap]()
  for k, thisheap in tables.mpairs(res):
    echo $k
    var newheap: ref MyHeap
    if not bread_to_areads.hasKey(k):
      new(newheap)
      newheap[] = binaryheap.newHeap[HeapItem](cmp_heap_items)
      bread_to_areads[k] = newheap
    else:
      newheap = bread_to_areads[k]
    for item in binaryheap.items(thisheap[]):
      if binaryheap.size(newheap[]) < settings.bestn:
        binaryheap.push(newheap[], item)
      else:
        discard binaryheap.pushpop(newheap[], item)
  #[
  For each b-read, we find the best contig map throgh the b->a->contig map.
  ]#
  let rawread_to_contigs_fn = settings.output
  let out_f = open(rawread_to_contigs_fn, fmWrite)
  defer: out_f.close()
  for bread in tables.keys(bread_to_areads):
    type MyArray = array[2, int]
    let ctg_score = tables.newTable[string, array[2, int]]()
    for s, rid in binaryheap.items(bread_to_areads[bread][]):
      if not rid_to_ctg.hasKey(rid): continue
      let ctgs = rid_to_ctg[rid]
      for ctg in ctgs:
        discard tables.mgetOrPut(ctg_score, ctg, [0,0])
        ctg_score[ctg][0] += -s
        ctg_score[ctg][1] += 1
    ########oid = rid_to_oid[int(bread)]
    type ctg_score_pair = tuple[ctg: string, score_count: array[2, int]]
    var ctg_scores: seq[ctg_score_pair] = sequtils.toSeq(tables.pairs(ctg_score))
    algorithm.sort(ctg_scores) do (x, y: ctg_score_pair) -> int:
      result = x.score_count[0]

    var rank = 0
    for cp in ctg_scores:
      var in_ctg = 0
      if rid_to_ctg.hasKey(bread):
        if rid_to_ctg[bread].find(cp.ctg) != -1:
          in_ctg = 1
      let score = cp.score_count[0]
      let count = cp.score_count[1]
      ######print(bread, oid, ctg, count, rank, score, in_ctg, file=out_f)
      let line = "$# $# $# $# $# $#" % [
          bread, cp.ctg, $count, $rank, $score, $in_ctg]
      out_f.write(line)
      rank += 1
  #[
  inputs = []
  for fn in file_list:
      inputs.append( (run_tr_stage1, db_fn, fn, min_len, bestn, rid_to_ctg, rid_to_phase) )
  bread_to_areads = {}
  for fn, res in exe_pool.imap(io.run_func, inputs):
      for k in res:
          bread_to_areads.setdefault(k, [])
          for item in res[k]:
              if len(bread_to_areads[k]) < bestn:
                  heappush( bread_to_areads[k], item )
              else:
                  heappushpop( bread_to_areads[k], item )

  #rid_to_oid = open(os.path.join(rawread_dir, 'dump_rawread_ids', 'rawread_ids')).read().split('\n')
  ]#
proc try_run_track_reads(s: Settings) =
  try:
    run_track_reads(s)
  except Exception:
    raise
  discard
proc track_reads(s: Settings) =
  echo "HI"
  try_run_track_reads(s)
proc main() =
  var settings: Settings
  track_reads(settings)

proc test() =
  block:
    let data_get_rid_to_ctg = """

  000000026 000000026 ref1/27/0_39949 000000F
  000000026 000000026 ref1/27/0_39949 000000C
  000000278 000000287 ref2/128/0_42163 000000F_02
  """
    let res = stream_get_rid_to_ctg(streams.newStringStream(
        data_get_rid_to_ctg))
    assert len(res) == 2, $len(res)
    assert $res["000000026"] == "@[000000F, 000000C]", $res["000000026"]
    assert $res["000000287"] == "@[000000F_02]", $res["000000287"]

if isMainModule:
  test()
  main()
