import threadpool

proc foo() =
    echo "foo"

proc main(min_cov=6, min_cov_aln=10, max_cov_aln=0, min_len_aln=0, min_n_read=10, max_n_read=500,
          trim=false, output_full=false, output_multi=false,
          min_idt="0.70", edge_tolerance=1000, trim_size=50,
          n_core=24): int =
    if n_core == 0:
        foo()
    else:
        spawnX foo()
    0

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
