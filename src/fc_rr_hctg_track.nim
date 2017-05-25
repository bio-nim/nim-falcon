# vim: sw=2 ts=2 sts=2 tw=80 et:
#import nimprof
from falcon/sys import log
from falcon/rr_hctg_track import nil

proc rr_hctg_track1(
    test=false,
    phased_read_file="./3-unzip/all_phased_reads",
    read_to_contig_map="./4-quiver/read_maps/read_to_contig_map",
    rawread_ids="./2-asm-falcon/read_maps/dump_rawread_ids/rawread_ids",
    bestn=40,
    n_core=0,
    min_len=2500,
    # All the rest are ignored....
    output="./2-asm-falcon/read_maps/dump_rawread_ids/rawread_to_contigs",
    stream=false, debug=false, silent=false,
  ) =
  if test:
    log("no tests")
    system.quit(system.QuitSuccess)
  let db_fn = "0-rawreads/raw_reads.db"
  let file_list = rr_hctg_track.glob_raw_reads_las()
  let settings = rr_hctg_track.Settings(
    n_core: n_core,
    file_list: file_list,
    db_fn: db_fn,
    phased_read_file: phased_read_file,
    read_to_contig_map: read_to_contig_map,
    rawread_ids: rawread_ids,
    output: output,
    min_len: min_len,
    stream: stream,
    debug: debug,
    bestn: bestn,
  )
  rr_hctg_track.run_track_reads(settings)
  echo "BYE"


when isMainModule:
  import "cligen/cligen"
  # cligen supports --foo-bar and --foo_bar styles automatically.
  dispatch(rr_hctg_track1)
