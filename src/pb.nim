# vim: sw=2 ts=2 sts=2 tw=0 et:
#from fc_consensus import nil
from fc_rr_hctg_track import nil
from falcon/rr_hctg_track import nil
from falcon/plasmid import nil

proc dataset(extras: seq[string]) =
  echo "pb dataset"

when isMainModule:
  import cligen
  dispatch(
        dataset, help = {},
  )
#[
  dispatchMulti(
        [dataset, help = {}],
        [plasmid.main, cmdName="plasmid"],
        # [kmers, short = {"int_dummy": 'd'}, help = {}],
        [fc_rr_hctg_track.rr_hctg_track1, cmdName="rr_hctg_track"],
        [rr_hctg_track.run_stage2, cmdName="rr_hctg_track2"],
        [fc_consensus.main, cmdName="consensus",
          short={"": '\0'},
          help={
    "min-idt": "minimum identity of the alignments used for correction (32-bit float)",
    "edge-tolerance": "for trimming, the there is unaligned edge leng > edge_tolerance, ignore the read",
    "trim": "trim the input sequence with k-mer spare dynamic programming to find the mapped range",
    "trim-size": "the size for triming both ends from initial sparse aligned region",
    "allow-external_mapping": "if provided, externally determined mapping coordinates will be used for error correction",
    "output-multi": "output multiple correct regions",
    "output-full": "output uncorrected regions too",
    "min-n-read": "1 + minimum number of reads used in generating the consensus; a seed read with fewer alignments will be completely ignored",
    "max-n-read": "1 + maximum number of reads used in generating the consensus",
    "n-core": "number of processes used for generating consensus (not sure this limit works yet); 0 for main process only",
    "min-cov": "minimum coverage to break the consensus",
    "min-cov-aln": "minimum coverage of alignment data; a seed read with less than MIN_COV_ALN average depth of coverage will be completely ignored",
    "max-cov-aln": "maximum coverage of alignment data; a seed read with more than MAX_COV_ALN average depth of coverage of the longest alignments will be capped, excess shorter alignments will be ignored",
    "min-len-aln": "minimum length of a sequence in an alignment to be used in consensus; any shorter sequence will be completely ignored",
        }],
    )
]#
