# vim: sw=2 ts=2 sts=2 tw=80 et:
from strutils import nil

proc padLeft*(s: string, width: Natural): string =
  # https://forum.nim-lang.org/t/1235/2
  if width > len(s):
    result = strutils.spaces(width - len(s)) & s
  else:
    result = s
