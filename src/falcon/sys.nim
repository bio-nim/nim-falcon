# vim: sw=2 ts=2 sts=2 tw=80 et:
const version* = "0.0.0"

from functional import padLeft

from cpuinfo import nil
#from locks import nil
from os import nil
from posix import nil
from strutils import nil
from threadpool import nil
from times import nil

proc log*[Ty](x: varargs[Ty, `$`]): auto =
  let start {.global.} = times.epochTime()
  let msg = strutils.join(x, "") # TODO: add spaces here, not at call-site
  let sofar = (times.epochTime() - start)
  let formatted = padLeft(strutils.formatFloat(sofar, format=strutils.ffDecimal, precision=1), 7)
  writeLine(stderr, formatted, "s ", msg)

# Show version at start-up, for now.
log("version=", $version, ", nim-version=", system.NimVersion)

proc abspath*(fn: string): auto =
  if not strutils.startswith(fn, '/'):
    return os.getCurrentDir() & '/' & fn
  else:
    return fn

proc filelength*(fn: string): auto =
  var
    st: posix.Stat
  let ret = posix.stat(fn.cstring, st)
  if 0.cint != ret:
    let code = os.osLastError()
    let msg = "Cannot stat('" & abspath(fn) & "')"
    os.raiseOSError(code, msg)
  assert ret == 0.cint
  return st.st_size

#[
var the_wait_is_over: locks.Cond
var the_wait_is_over_lock: locks.Lock

proc waiter() =
  log("Wait...")
  locks.wait(the_wait_is_over, the_wait_is_over_lock)
  log("Stopped waiting.")

proc setMaxPoolSize*(nmax: int) =
  let currMax = threadpool.MaxThreadPoolSize
  threadpool.setMaxPoolSize(nmax)
  locks.initLock(the_wait_is_over_lock)
  locks.initCond(the_wait_is_over)
  defer: locks.deinitLock(the_wait_is_over_lock)
  defer: locks.deinitCond(the_wait_is_over)
  log("Creating ", currMax, " dummies...")
  var count = 0
  for i in 0..<24: #currMax:
    if not threadpool.preferSpawn():
      break
    threadpool.spawn waiter()
    count.inc
  log("Created ", count, " waiters. Signaling and syncing...")
  while count > 0:
    locks.signal(the_wait_is_over)
    count.dec
  threadpool.sync()
  log("threadpool.MaxThreadPoolSize was:", currMax, " now:", threadpool.MaxThreadPoolSize, ", should be:", nmax)
]#

proc dummy() {.thread.} =
  os.sleep(500)

proc setMaxPoolSize*(nmax: int) =
  threadpool.setMaxPoolSize(nmax)
  let ncpus = cpuinfo.countProcessors()
  log("Creating ", ncpus, " dummies...")
  for i in 0..<ncpus:
    threadpool.spawn dummy()
  threadpool.sync()
  log(" Finished all dummies. The threadpool size should now be ", nmax)
