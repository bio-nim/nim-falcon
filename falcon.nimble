# Package

version       = "0.0.0"
author        = "Christopher Dunn"
description   = "Nim code for running Falcon genome assembler."
license       = "BSD 3 Clear"

# Dependencies

requires "nim >= 0.19.9", "daligner", "binaryheap", "cligen", "htslib", "msgpack4nim"

srcDir = "./src"
bin = @["falcon/example", "fc_rr_hctg_track2", "fc_rr_hctg_track", "fc_consensus"]

#task test, "Test daligner wrapper":
#    withDir("tests"):
#        exec("make")
