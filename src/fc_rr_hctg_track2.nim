# vim: sw=2 ts=2 sts=2 tw=80 et:
#import nimprof
from falcon/rr_hctg_track import nil


when isMainModule:
  import "cligen"
  # cligen supports --foo-bar and --foo_bar styles automatically.
  dispatch(rr_hctg_track.run_stage2)
