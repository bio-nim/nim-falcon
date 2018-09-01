# vim: sw=2 ts=2 sts=2 tw=0 et:
#import nimprof
from json import `[]`
from os import nil
from ospaths import nil
from posix import nil
from strformat import fmt
from strtabs import nil
from tables import nil
const DONE = "done"
const STATUS = "status"
var TIMEOUT = 1 # global and default
const DESCRIPTION = """Given a JSON description, call a python-function.
"""
const EPILOG = """
The JSON looks like this:
{
    "inputs": {"input-name": "filename"},
    "outputs": {"output-name": "output-filename (relative)"},
    "bash_template_fn": "template.sh",
    "parameters": {}
}

This program will run on the work host, and it will do several things:
    - Run in CWD.
    - Verify that inputs are available. (Wait til timeout if not.)
    - Possibly, cd to tmpdir and create symlinks from inputs.
    - Run the python-function.
      - Its module must be available (e.g. in PYTHONPATH).
      - Pass a kwd-dict of the union of inputs/outputs/parameters.
      - Ignore return-value. Expect exceptions.
    - Possibly, mv outputs from tmpdir to workdir.
    - Write exit-code into STATUS.
    - Touch DONE on success.
"""
discard """
(Someday, we might also support runnable Python modules, or even executables via execvp().)

Note: qsub will *not* run this directly. There is a higher layer.
"""

template with_cd*(d: string, body: untyped) =
  block:
    var old_dir = os.getCurrentDir()
    os.setCurrentDir(d)
    echo "<-", d
    body
    echo "->", old_dir
    os.setCurrentDir(old_dir)
proc wait_for(fn: string) =
    var timeout = TIMEOUT
    #LOG.debug('Checking existence of {!r} with timeout={}'.format(fn, timeout))
    var dirname = ospaths.parentDir(fn)
    if os.existsFile(dirname):
      if false: #not os.access(dirname, os.X_OK):
        raise newException(ValueError, "Cannot x into dir {dirname}".fmt)
    while not os.existsFile(fn):
      if timeout > 0:
        os.sleep(1000)
        timeout -= 1
      else:
        raise newException(ValueError, "Timed out waiting for {fn}".fmt)
    #assert os.access(fn, os.R_OK), '{!r} not readable'.format(fn)
discard """
class Attrs(object):
    #This facilitates substitution of values in string.
    def __str__(self):
        # For this, all values must be strings.
        return ' '.join(f for f in self.kwds.values())
    def __getattr__(self, name):
        # For this, values can be string, int, float, etc.
        return str(self.kwds[name])
    def __init__(self, **kwds):
        self.kwds = kwds
proc value_quoted(kvs):
    return {k:quoteShell(v) for k,v in kvs.items()}
"""
proc run_bash(bash_template: string, myinputs, myoutputs, parameters: tables.OrderedTable[string, string]) =
  discard """
    # Set substitution dict
    var_dict = dict()
    #var_dict.update(parameters)
    #var_dict.update(myinputs) # for now
    #var_dict.update(myoutputs) # for now
    valid_parameters = {k:v for k,v in parameters.iteritems() if not k.startswith('_')}
    assert 'input' not in parameters
    assert 'output' not in parameters
    # input/output/params are the main values substituted in the subset of
    # snakemake which we support.
    var_dict['input'] = Attrs(**value_quoted(myinputs))
    var_dict['output'] = Attrs(**value_quoted(myoutputs))
    var_dict['params'] = Attrs(**valid_parameters)
  """
  # Like snakemake, we use bash "strict mode", but we add -vx.
  # http://redsymbol.net/articles/unofficial-bash-strict-mode/
  var prefix = """
IFS=$'\n\t'
set -vxeuo pipefail
hostname
pwd
date
"""
  var postfix = """
date
"""
  # Substitute
  var bash_content = ""
  try:
      # Before subst, we need to add $ before {}s.
      # ...
      # Now,
      bash_content = prefix & postfix #& bash_template.format(**var_dict) & postfix
  except:
      var msg = """\
Failed to substitute var_dict
{}
into bash script:
{}
Possibly you forgot to use "input.foo" "output.bar" "params.fubar" etc. in your script?
"""#.format(var_dict, bash_template)
      #LOG.error(msg)
      raise
  discard """
    # Write user_script.sh
    bash_fn = 'user_script.sh'
    with open(bash_fn, 'w') as ofs:
        ofs.write(bash_content)
    cmd = '/bin/bash {}'.format(bash_fn)
    util.system(cmd)
  """
proc run_cfg_in_tmpdir(cfg: json.JsonNode, tmpdir: string) =
  echo "Running?"
  # Accept 'inputs', 'outputs', 'parameters' in cfg.
  # Require 'bash_template_fn' in cfg; substitute and use it.
  var
    myinputs = json.to(cfg["inputs"], tables.OrderedTable[string,string])
    myoutputs = json.to(cfg["outputs"], tables.OrderedTable[string,string])
    parameters = json.to(cfg["inputs"], tables.OrderedTable[string,string])
    bash_template_fn = json.getStr(cfg["bash_template_fn"])
    myrundir = ""
  for fn in tables.values(myinputs):
      echo "fn", fn
      wait_for(fn)
  wait_for(bash_template_fn)
  var bash_template = system.readFile(bash_template_fn)
  echo bash_template
  var finaloutdir = os.getCurrentDir()
  if tmpdir != "":
    discard """
      import getpass
      user = getpass.getuser()
      pid = os.getpid()
      myrundir = '{tmpdir}/{user}/pypetmp/{finaloutdir}'.format(**locals())
      util.rmdirs(myrundir)
      util.mkdirs(myrundir)
      # TODO(CD): Copy inputs w/ flock.
    """
  else:
    myrundir = finaloutdir
  with_cd(myrundir):
    # TODO(CD): Write a script in wdir even when running in tmpdir.
    run_bash(bash_template, myinputs, myoutputs, parameters)
  discard """
  if tmpdir:
      cmd = 'rsync -av {}/ {}; rm -rf {}'.format(myrundir, finaloutdir, myrundir)
      util.system(cmd)
  """
  for fn in tables.values(myoutputs):
      wait_for(fn)
proc normPath(p: string): string =
  if p == "":
    return "."
  return p
proc main(timeout: int=TIMEOUT, tmpdir: string=nil, json_fn: seq[string]) =
  echo "hi"
  wait_for(json_fn[0])
  TIMEOUT = timeout
  #LOG.debug('Loading JSON from {!r}'.format(json_fn[0]))
  var cfg = json.parseFile(json_fn[0])
  echo json.pretty(cfg)
  #LOG.debug(pprint.pformat(cfg))
  var rundir = normPath(ospaths.parentDir(json_fn[0]))
  echo "rundir:", rundir
  with_cd(rundir):
    echo "cwd:", os.getCurrentDir()
    run_cfg_in_tmpdir(cfg, tmpdir)
when isMainModule:
  import "cligen/cligen"
  cligen.dispatch(main, short={}, help={
    "timeout": "How many seconds to wait for input files (and JSON) to exist.",
    "tmpdir": "Root directory to run in. (Sub-dir name will be based on CWD.)",
  })
