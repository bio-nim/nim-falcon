# vim: sw=2 ts=2 sts=2 tw=80 et:
#import nimprof

from falcon/consensus/common import log
from falcon/consensus/falcon import nil
from falcon/consensus/DW_banded import nil
from falcon/consensus/kmer_lookup_c import nil
#from falcon/consensus/poo import nil

import os
import algorithm
import threadpool
import sequtils
import sets
from strutils import `strip`, `split`
from strutils import `%`, `formatFloat`, `ffdecimal`

# This has been problematic.
# But the main reason to avoid Regex is described here:
#   https://nim-lang.org/blog/2015/10/27/version-0120-released.html
#   "could not import: pcre_free_study"
when defined(use_pcre):
  import nre
  var good_regions {.threadvar.}: Regex

proc get_longest_reads(seqs: seq[string], max_n_read, max_cov_aln: int): seq[string] =
    var longest_n_reads = max_n_read
    if max_cov_aln > 0:
        longest_n_reads = 1
        var seed_len = len(seqs[0])
        var read_cov = 0
        for seq in seqs[1 .. seqs.high]:
            if read_cov div seed_len > max_cov_aln:
                break
            longest_n_reads += 1
            read_cov += len(seq)

        longest_n_reads = min(longest_n_reads, max_n_read)

    #assert longest_n_reads <= len(seqs), "{} <= {}".format(longest_n_reads, len(seqs))
    longest_n_reads = min(longest_n_reads, len(seqs)) # unnec in python
    return seqs[0 .. longest_n_reads-1]
proc get_longest_sorted_reads(seqs: seq[string], max_n_read, max_cov_aln: int): seq[string] =
  # allows us to avoid a redundant sort
  # in get_consensus_trimmed()
  var sorted_seqs: seq[string] = @[]
  sorted_seqs.insert(seqs[1 .. seqs.high], 0)
  #seqs[1 .. seqs.high].sort(proc (x,y: string): int = cmp(y.len, x.len))
  sorted_seqs.sort(proc (x,y: string): int = cmp(y.len, x.len))
  sorted_seqs.insert([seqs[0]], 0)
  #log("sortlen", $len(sorted_seqs))
  return get_longest_reads(sorted_seqs, max_n_read, max_cov_aln)
type
  Config = tuple
    output_multi: bool
    min_cov: int
    K: int
    max_n_read: int
    min_idt: float32
    edge_tolerance: int
    trim_size: int
    allow_external_mapping: bool
    min_cov_aln: int
    max_cov_aln: int
iterator get_seq_data(config: Config, min_n_read, min_len_aln: int): auto =
  #let min_idt = 0.70
  #let min_cov = 1
  #let max_n_read 20000
  const max_len = 100000
  var seed_len = 0
  var seed_id: string = ""
  var seqs: seq[string] = @[]
  var read_ids: HashSet[string]
  var read_cov: int
  init(read_ids)
  var inp = stdin
  #const fn = "/Users/cdunn2001/repo/gh/nim-consensus/c/t/foo.long"
  #var inp = open(fn)
  for fullline in inp.lines:
    var split_line = fullline.strip.split()
    if len(split_line) != 2: continue
    var read_id = split_line[0]
    var sequ = split_line[1]
    if len(sequ) > max_len:
        sequ = sequ[0 .. max_len-1]
    case read_id
    of "+":
      if len(seqs) >= min_n_read and read_cov div seed_len >= config.min_cov_aln:
          #log("len(seqs) ", $len(seqs))
          seqs = get_longest_sorted_reads(seqs, config.max_n_read, config.max_cov_aln)
          yield (seqs, seed_id, config)
      seqs = @[]
      read_ids.init
      seed_id = ""
      read_cov = 0
    of "-":
      break
    of "*":
      seqs = @[]
      read_ids.init
      seed_id = ""
      read_cov = 0
    else:
      if len(sequ) >= min_len_aln:
          if len(seqs) == 0:
              seqs.add(sequ) #the "seed"
              seed_len = len(sequ)
              seed_id = read_id
          if not (read_id in read_ids): #avoidng using the same read twice. seed is used again here by design
              seqs.add(sequ) #the "seed"
              read_ids.incl(read_id)
              read_cov += len(sequ)
proc copy_seq_ptrs(cseqs: var cStringArray, seqs: seq[string]) =
  cseqs = allocCStringArray(seqs)

type ConsensusArgs = tuple
  inseqs: seq[string]
  seed_id: string
  config: Config
type ConsensusResult = tuple
  consensus: string
  seed_id: string
proc get_con(args: ConsensusArgs): ConsensusResult =
    var config = args.config
    var seqs = args.inseqs
    let n_seq = len(seqs)
    if len(seqs) > config.max_n_read:
        seqs = get_longest_sorted_reads(seqs, config.max_n_read, config.max_cov_aln)
    #poo.poo()
    #log("About to generate_con ", $len(seqs), " ", $n_seq)
    #log("pseq:", $(cast[ByteAddress]((addr(seqs[0])))), " eg ", $len(seqs[0]))
    var cseqs: cStringArray
    copy_seq_ptrs(cseqs, seqs)
    var consensus_data = falcon.generate_consensus(cseqs, n_seq, config.min_cov, config.K, config.min_idt)
    deallocCStringArray(cseqs)
    #log("cr:", repr(consensus_data))
    var consensus = consensus_data.sequence
    #echo "cs:", addr(consensus_data.sequence), " ", addr consensus, " ", addr consensus_data.sequence[0], " ", addr consensus[0]
    #eff_cov = consensus_data_ptr.eff_cov[:len(consensus)]
    return (consensus, args.seed_id)
discard """
proc get_consensus_without_trim(args: ConsensusArgs): ConsensusResult =
    var config = args.config
    var seqs = args.inseqs
    if len(seqs) > config.max_n_read:
        seqs = get_longest_sorted_reads(seqs, config.max_n_read, config.max_cov_aln)
    #seqs_ptr = seqs # copy
    #for i in countup(0, len(seqs)-1):
    #  log("i len(seq) %d %d"%(i, len(seq)))
    var cseqs: cStringArray
    let n_seq: cuint = cuint(len(seqs))
    copy_seq_ptrs(cseqs, seqs)
    #echo cseqs

    var all_seqs_mapped = false

    if config.allow_external_mapping:
      all_seqs_mapped = true
      for seq in seqs:
      if not seq.is_mapped:
        all_seqs_mapped = false
        break

    if not all_seqs_mapped:
      #foo.poo()
      log("Internally mapping the sequences.")
      log("About to generate_consensus ", $len(seqs), " ", $n_seq)
      consensus_data_ptr = common.generate_consensus(
          cseqs,
          n_seq, cuint(config.min_cov), cuint(config.K), cdouble(config.min_idt))
    else:
      log("Using external mapping coordinates from input.")
      var aln_ranges = newSeq[aln_range](len(seqs))
      for i in [0 ..< len(seqs)]:
        let seq = seqs[i]
        aln_ranges[i] = aln_range(
            s1: seq.qstart, e1: seq.qend,
            s2: seq.tstart, e2: seq.tend,
            score: (seq.qend - seq.qstart))
      log("About to generate_consensus ", $len(seqs), " ", $n_seq)
      consensus_data_ptr = common.generate_consensus_from_mapping(
          cseqs,
          aln_ranges,
          n_seq, cuint(config.min_cov), cuint(config.K), cdouble(config.min_idt))

    deallocCStringArray(cseqs)
    var consensus = $consensus_data_ptr.sequence # expect a copy
    #eff_cov = consensus_data_ptr.eff_cov[:len(consensus)]
    common.free_consensus_data(consensus_data_ptr)
    return (consensus, args.seed_id)
  """

when defined(use_pcre):
  proc findall_patt(consensus: string, patt: Regex): seq[string] =
    result = findall(consensus, patt)
    #echo consensus[0], " ", len(consensus), " ", consensus[^1], " ", len(result)

iterator findall_good_regions(consensus: string, outoo: var string): int =
  let n = len(consensus)
  var cbeg = 0
  var cend = 0
  while cend < n:
    let c = consensus[cend]
    if c == 'A' or c == 'C' or c == 'G' or c == 'T':
      inc(cend)
      continue
    if cend != cbeg:
      outoo.setLen(cend-cbeg)
      outoo[0..^1] = consensus[cbeg..<cend]
      yield (cend-cbeg)
    inc(cend)
    cbeg = cend
  if cend != cbeg:
    outoo.setLen(cend-cbeg)
    outoo[0..^1] = consensus[cbeg..<cend]
    yield (cend-cbeg)

proc format_seq(sequ: string, col: int): string =
  result = newString(len(sequ) + len(sequ) div col + 1)
  var bo = 0
  var bn = 0
  while (bo+col) < len(sequ):
    result[bn ..< (bn+col)] = sequ[bo ..< (bo+col)]
    result[(bn+col)] = '\l'
    bo += col
    bn += col + 1
  var tail = len(sequ) - bo
  result[bn ..< (bn+tail)] = sequ[bo ..< (bo+tail)]
  result.setLen(bn+tail)
  #result[(bn+tail)] = '\l' # Python did not add final newline
proc write_seq(cns_seq: string, seed_id: string, seq_i: var int): bool =
            # Return false when done.
            if len(cns_seq) < 500:
              return true
            if seq_i >= 10:
              return false
            #print ">prolog/%s%01d/%d_%d" % (seed_id, seq_i, 0, len(cns_seq))
            echo ">prolog/", seed_id, seq_i, "/", 0, "_", len(cns_seq)
            echo format_seq(cns_seq, 80)
            seq_i += 1
            return true
proc process_consensus(cargs: ConsensusArgs) {.thread} =
    #discard GC_disable
    discard """
    if thread_msa_array == nil:
      log "Was nil"
      thread_msa_array = ""
    else:
      log "Was not nil"
    """
    #var (consensus, seed_id) = get_consensus_without_trim(cargs)
    var (consensus, seed_id) = get_con(cargs)
    let config = cargs.config
    #log("len(consensus)=", $len(consensus), " in seed ", seed_id)
    if len(consensus) < 500:
        return
    if false: # args.output_full:
        echo ">"&seed_id&"_f"
        echo consensus
        return
    when defined(use_pcre):
      var good_regions: Regex
      if good_regions.isNil:
        good_regions = re"[ACGT]+"
        #  #log("good_regions isNil!!!")
      #var cns = findall_patt(consensus, good_regions)
    if config.output_multi:
        var seq_i = 0
        when defined(use_pcre):
          for cns_seq in findall_patt(consensus, good_regions):
            let more = write_seq(cns_seq, seed_id, seq_i)
        else:
          var cns_seq = newStringOfCap(len(consensus))
          for _ in findall_good_regions(consensus, cns_seq):
            #log("$# $#\L $#" % [$len(cns_seq), $len(consensus), repr(cns_seq)]) #, repr(consensus)])
            let more = write_seq(cns_seq, seed_id, seq_i)
            if not more:
              break
    else:
        raiseAssert("--output-multi is required for now.")
    #    var cns = findall_good_regions(consensus)
    #    #cns.sort(key = lambda x: len(x))
    #    echo ">"&seed_id
    #    echo cns[cns.high-1]
discard """
proc simple(cargs: ConsensusArgs) =
  echo "hi in simple()"
  #var (consensus, seed_id) = get_consensus_without_trim(cargs)
  discard get_consensus_without_trim(cargs)
  """

proc waitForOpenThreadIndex(threads: var seq[ref Thread[ConsensusArgs]]): int =
  let n = len(threads)
  assert(n > 0)
  for i in 0..<n:
    if nil == threads[i]:
      return i
  while true:
    for i in 0..<n:
      if not running(threads[i][]):
        return i
    os.sleep(100)
proc main(min_cov=6, min_cov_aln=10, max_cov_aln=0, min_len_aln=0, min_n_read=10, max_n_read=500,
          trim=false, output_full=false, output_multi=false,
          min_idt="0.70", edge_tolerance=1000, trim_size=50, allow_external_mapping=false,
          n_core=24): int =
  doAssert(output_multi, "--output-multi is required for now.")
  doAssert(not trim, "--trim is required for now.")
  log("main(n_core=", $n_core, ")")
  if n_core > 0:
    #threadpool.setMaxPoolSize(n_core)
    #sync() # Let extra threads shutdown.
    discard
  let config: Config = (
    output_multi: output_multi,
    min_cov: min_cov,
    K: 8, # not cli
    max_n_read: max_n_read,
    min_idt: float32(strutils.parseFloat(min_idt)),
    edge_tolerance: edge_tolerance,
    trim_size: trim_size,
    allow_external_mapping: allow_external_mapping,
    min_cov_aln: min_cov_aln,
    max_cov_aln: max_cov_aln)
  #log("config=", config)
  var threads = newSeq[ref Thread[ConsensusArgs]](n_core)
  for q in get_seq_data(config, min_n_read, min_len_aln):
    var (seqs, seed_id, config_same) = q
    log("len(seqs)=", $len(seqs), ", seed_id=", seed_id)
    var cargs: ConsensusArgs = (inseqs: seqs, seed_id: seed_id, config: config)
    if n_core == 0:
      #common.benchmark "loop":
      process_consensus(cargs)
    else:
      #let thread = wait(threads)
      var rthread: ref Thread[ConsensusArgs]
      new(rthread)
      let i = waitForOpenThreadIndex(threads)
      #threads.add(rthread)
      threads[i] = rthread
      createThread(rthread[], process_consensus, cargs)
      #joinThread(rthread[])
    #  spawn process_consensus(cargs)
    #  #spawn os.sleep(1000)
    #  #spawn simple(cargs)
    #log("tot=$1 occ=$2, free=$3 b4" % [$getTotalMem(), $getOccupiedMem(), $getFreeMem()])
    #GC_fullCollect()
    #log("tot=$1 occ=$2, free=$3 now" % [$getTotalMem(), $getOccupiedMem(), $getFreeMem()])
  #sync()
  #log("Num threads still running:", $len(threads)) # not necessarily
  for rthread in threads:
    if nil != rthread:
      joinThread(rthread[])
      log("Finished a thread.")
  log("tot=$1 occ=$2, free=$3 b4" % [$getTotalMem(), $getOccupiedMem(), $getFreeMem()])
  GC_fullCollect()
  log("tot=$1 occ=$2, free=$3 now" % [$getTotalMem(), $getOccupiedMem(), $getFreeMem()])
  result = 0

when isMainModule:
  import "cligen/cligen"
  cligen.dispatch(main, short={}, help={
    "min_idt": "minimum identity of the alignments used for correction (32-bit float)",
    "edge_tolerance": "for trimming, the there is unaligned edge leng > edge_tolerance, ignore the read",
    "trim": "trim the input sequence with k-mer spare dynamic programming to find the mapped range",
    "trim_size": "the size for triming both ends from initial sparse aligned region",
    "allow_external_mapping": "if provided, externally determined mapping coordinates will be used for error correction",
    "output_multi": "output multiple correct regions",
    "output_full": "output uncorrected regions too",
    "min_n_read": "1 + minimum number of reads used in generating the consensus; a seed read with fewer alignments will be completely ignored",
    "max_n_read": "1 + maximum number of reads used in generating the consensus",
    "n_core": "number of processes used for generating consensus (not sure this limit works yet); 0 for main process only",
    "min_cov": "minimum coverage to break the consensus",
    "min_cov_aln": "minimum coverage of alignment data; a seed read with less than MIN_COV_ALN average depth of coverage will be completely ignored",
    "max_cov_aln": "maximum coverage of alignment data; a seed read with more than MAX_COV_ALN average depth of coverage of the longest alignments will be capped, excess shorter alignments will be ignored",
    "min_len_aln": "minimum length of a sequence in an alignment to be used in consensus; any shorter sequence will be completely ignored",
  })
