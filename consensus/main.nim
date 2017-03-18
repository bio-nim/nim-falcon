# vim: sw=2 ts=2 sts=2 tw=80 et:
#import nimprof

from common import log
from falcon import nil
from DW_banded import nil
from kmer_lookup_c import nil
from poo import nil

import os
import algorithm
import threadpool
import nre
import sequtils
import sets
from strutils import `strip`, `split`
from strutils import `%`, `formatFloat`, `ffdecimal`


let good_regions = re"[ACGT]+"
#var thread_msa_array {.threadvar.}: string

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
    min_cov: int
    K: int
    max_n_read: int
    min_idt: float32
    edge_tolerance: int
    trim_size: int
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
    var line = fullline.strip.split()
    if len(line) != 2: continue
    var read_id = line[0]
    var sequ = line[1]
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
  good_regions: Regex
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
    #foo.poo()
    log("About to generate_consensus ", $len(seqs), " ", $n_seq)
    #echo cseqs
    var consensus_data_ptr = common.generate_consensus(cseqs, n_seq, cuint(config.min_cov), cuint(config.K), cdouble(config.min_idt))
    deallocCStringArray(cseqs)
    var consensus = $consensus_data_ptr.sequence # expect a copy
    #eff_cov = consensus_data_ptr.eff_cov[:len(consensus)]
    common.free_consensus_data(consensus_data_ptr)
    return (consensus, args.seed_id)
  """
proc findall_patt(consensus: string, patt: Regex): seq[string] =
  result = findall(consensus, patt)
  #echo consensus[0], " ", len(consensus), " ", consensus[^1], " ", len(result)
proc format_seq(sequ: string, col: int): string =
  result = newString(len(sequ) + len(sequ) div col + 1)
  var bo = 0
  var bn = 0
  while (bo+col) < len(sequ):
    result[bn .. <(bn+col)] = sequ[bo .. <(bo+col)]
    result[(bn+col)] = '\l'
    bo += col
    bn += col + 1
  var tail = len(sequ) - bo
  result[bn .. <(bn+tail)] = sequ[bo .. <(bo+tail)]
  result.setLen(bn+tail)
  #result[(bn+tail)] = '\l' # Python did not add final newline
proc process_consensus(cargs: ConsensusArgs) {.thread} =
    discard """
    if thread_msa_array == nil:
      log "Was nil"
      thread_msa_array = ""
    else:
      log "Was not nil"
    """
    #var (consensus, seed_id) = get_consensus_without_trim(cargs)
    var (consensus, seed_id) = get_con(cargs)
    #log("len(consensus)=", $len(consensus), " in seed ", seed_id)
    if len(consensus) < 500:
        return
    if false: # args.output_full:
        echo ">"&seed_id&"_f"
        echo consensus
        return
    var cns = findall_patt(consensus, cargs.good_regions)
    #var cns = findall_patt(consensus, good_regions)
    if len(cns) == 0:
        return
    if true: #args.output_multi:
        var seq_i = 0
        for cns_seq in cns:
            #log("%d %d" %(len(cns_seq), len(consensus)))
            if len(cns_seq) < 500:
                return
            if seq_i >= 10:
                break
            #print ">prolog/%s%01d/%d_%d" % (seed_id, seq_i, 0, len(cns_seq))
            echo ">prolog/", seed_id, seq_i, "/", 0, "_", len(cns_seq)
            echo format_seq(cns_seq, 80)
            seq_i += 1
    else:
        #cns.sort(key = lambda x: len(x))
        echo ">"&seed_id
        echo cns[cns.high-1]
discard """
proc simple(cargs: ConsensusArgs) =
  echo "hi in simple()"
  #var (consensus, seed_id) = get_consensus_without_trim(cargs)
  discard get_consensus_without_trim(cargs)
  """

proc main(min_cov=6, min_cov_aln=10, max_cov_aln=0, min_len_aln=0, min_n_read=10, max_n_read=500,
          trim=false, output_full=false, output_multi=false,
          min_idt="0.70", edge_tolerance=1000, trim_size=50,
          n_core=24): int =
  log("main(n_core=", $n_core, ")")
  if n_core > 0:
    #threadpool.setMaxPoolSize(n_core)
    #sync() # Let extra threads shutdown.
    discard
  let config: Config = (
    min_cov: min_cov,
    K: 8, # not cli
    max_n_read: max_n_read,
    min_idt: float32(strutils.parseFloat(min_idt)),
    edge_tolerance: edge_tolerance,
    trim_size: trim_size,
    min_cov_aln: min_cov_aln,
    max_cov_aln: max_cov_aln)
  for q in get_seq_data(config, min_n_read, min_len_aln):
    var (seqs, seed_id, config_same) = q
    #log("len(seqs)=", $(len(seqs), ", seed_id=", seed_id, "config=", config))
    var cargs: ConsensusArgs = (inseqs: seqs, seed_id: seed_id, config: config, good_regions: good_regions)
    if n_core == 0:
      #common.benchmark "loop":
        process_consensus(cargs)
    else:
      spawnX process_consensus(cargs)
      #spawn os.sleep(1000)
      #spawn simple(cargs)
  sync()
  result = 0

when isMainModule:
  import "cligen/cligen"
  cligen.dispatch(main, short={}, help={
    "min_idt": "minimum identity of the alignments used for correction (32-bit float)",
    "edge_tolerance": "for trimming, the there is unaligned edge leng > edge_tolerance, ignore the read",
    "trim": "trim the input sequence with k-mer spare dynamic programming to find the mapped range",
    "trim_size": "the size for triming both ends from initial sparse aligned region",
    "output_multi": "output multi correct regions",
    "output_full": "output uncorrected regions too",
    "min_n_read": "1 + minimum number of reads used in generating the consensus; a seed read with fewer alignments will be completely ignored",
    "max_n_read": "1 + maximum number of reads used in generating the consensus",
    "n_core": "number of processes used for generating consensus (not sure this limit works yet); 0 for main process only",
    "min_cov": "minimum coverage to break the consensus",
    "min_cov_aln": "minimum coverage of alignment data; a seed read with less than MIN_COV_ALN average depth of coverage will be completely ignored",
    "max_cov_aln": "maximum coverage of alignment data; a seed read with more than MAX_COV_ALN average depth of coverage of the longest alignments will be capped, excess shorter alignments will be ignored",
    "min_len_aln": "minimum length of a sequence in an alignment to be used in consensus; any shorter sequence will be completely ignored",
  })
