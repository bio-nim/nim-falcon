# vim: sw=2 ts=2 sts=2 tw=80 et:
{.passC: "-g -Wall -I.".}
{.compile: "DW_banded.c".}
{.compile: "poo.c".}
{.compile: "falcon.c".}
{.compile: "kmer_lookup.c".}
import os
import algorithm
import threadpool
import nre
import sequtils
import sets
import strutils

from foo import nil
from common import nil

let good_regions = re"[ACGT]+"
#var thread_msa_array {.threadvar.}: string

proc log(msgs: varargs[string]) =
  #for s in msgs:
  #  write(stderr, s)
  #write(stderr, '\l')
  return
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
  log("sortlen", $len(sorted_seqs))
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
  const max_len = 10000
  var seed_len = 0
  var seed_id: string = ""
  var seqs: seq[string] = @[]
  var read_ids: HashSet[string]
  var read_cov: int
  init(read_ids)
  for line in stdin.lines:
    var l = line.strip.split()
    if len(l) != 2: continue
    var read_id = l[0]
    var sequ = l[1]
    if len(sequ) > max_len:
        sequ = sequ[0 .. max_len-1]
    case read_id
    of "+":
      if len(seqs) >= min_n_read and read_cov div seed_len >= config.min_cov_aln:
          log("len(seqs) ", $len(seqs))
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
proc process_consensus(cargs: ConsensusArgs) = #{.thread} =
    discard """
    if thread_msa_array == nil:
      log "Was nil"
      thread_msa_array = ""
    else:
      log "Was not nil"
    """
    var (consensus, seed_id) = get_consensus_without_trim(cargs)
    log($len(consensus), " in seed ", seed_id)
    if len(consensus) < 500:
        return
    if false: # args.output_full:
        echo ">"&seed_id&"_f"
        echo consensus
        return
    var cns = findall_patt(consensus, good_regions)
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
proc simple(cargs: ConsensusArgs) =
  echo "hi"
  #var (consensus, seed_id) = get_consensus_without_trim(cargs)
  discard get_consensus_without_trim(cargs)
proc main() =
  log("Setting MaxPoolSize")
  log("hi")
  let config: Config = (
    min_cov: 1, # default 6
    K: 8, # not cli
    max_n_read: 500,
    min_idt: 0.70f32, # default also
    edge_tolerance: 1000,
    trim_size: 50,
    min_cov_aln: 10,
    max_cov_aln: 0)
  let min_n_read = 10
  let min_len_aln = 0
  for q in get_seq_data(config, min_n_read, min_len_aln):
    var (seqs, seed_id, config_same) = q
    log($(len(seqs), seed_id, config))
    var cargs: ConsensusArgs = (inseqs: seqs, seed_id: seed_id, config: config)
    spawn process_consensus(cargs)
    #process_consensus(cargs)
    #spawn os.sleep(1000)
    #spawn simple(cargs)
  sync()
when isMainModule:
  threadpool.setMaxPoolSize(1) #n_cores)
  #threadpool.Xsetup()
  main()
