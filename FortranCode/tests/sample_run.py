from __future__ import print_function
from os.path import join, exists
import shutil
import re

def copyfiles(inpath, outpath):
  """ Copies input files. """
  from shutil import copyfile

  for pattern in ["calculation", "approximation"]:
    infilename = join(inpath, "{0}_input.dat".format(pattern))
    outfilename = join(outpath, "{0}_input.dat".format(pattern))
    if exists(infilename):
      copyfile(infilename, outfilename)
  inrestartfile = join(inpath, "restart_Input.inf")
  outrestartfile = join(outpath, "restart.inf")
  if exists(inrestartfile):
    copyfile(inrestartfile, outrestartfile)

def read_procstat(path):
  result = []
  with open(path, 'r') as file: 
    for line in file:
      data = line.split()
      if data[0] == 'Overall':
        pass # skip header. We do this rather than file.readline
             # because the header might be missing
             # (i.e. resume with no preexisting output)
      elif data[0] == 'configuration':
        result.append([float(u) for u in data[1:]])
      else: result[-1].extend([float(u) for u in data])
  return result

def loadtxt(path):
  with open(path, 'r') as file:
    file.readline() 
    return [ [float(u) for u in line.split()] for line in file ] 


def numberofevents(path,filename):
  import re

  setupfile = join(path, filename)
  Nevents = 0
  with open(setupfile) as file:
    for line in file.readlines():
      match = re.match("^ Events occurred",line)
      if match:
        eventline = line.split()
        Nevents = int(eventline[2])
  return Nevents

def doubleWallTime(directory, filename):
  fullfile = join(directory, filename)
  tempfilename = fullfile + 'new'
  backupfilename = fullfile + 'backup'
  with open(fullfile) as oldfile, open(tempfilename,'w') as tempfile:
      for line in oldfile:
          matchObj = re.search(r'^(wall_time\s*)(\d*)', line, re.M)
          if(matchObj):
              tempfile.write('{0}{1}\n'.format(matchObj.group(1), 2*int(matchObj.group(2))))
          else:
              tempfile.write(line)
  shutil.move(fullfile,backupfilename)
  shutil.move(tempfilename,fullfile)

def test(inpath, execpath, testpath, resumetest=False, Neventstarget=0, execrunner=None):
  """ Runs simulation and gets output. """
  from os.path import join, exists
  from shutil import rmtree
  from os import makedirs, getcwd
  from subprocess import Popen, PIPE
  from sys import exit

  directory = join(getcwd(),'mesoapprox_testoutput',testpath)
  exectorun = get_exectorun(execpath, execrunner)

  if exists(directory):
    shutil.rmtree(directory)
  if not exists(directory):
    makedirs(directory)
  try:
    datapath = join(inpath, join('data', testpath))
    copyfiles(datapath, directory)
    process = Popen(exectorun, cwd=directory, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

    print("Stdout from this test was:")
    print(stdout)
    print("Stderr from this test was:")
    print(stderr)
    
    fcalcin = open(join(directory,'calculation_input.dat'), 'r')
    # use readline() to read the first line which contains the approximation name
    line = fcalcin.readline()
    
    s = ''
    testdata = loadtxt(s.join((datapath, '\\', line, '_Fortran_Theta_vs_Mu.txt')))
    outdata = loadtxt(s.join((directory, '\\', line, '_Fortran_Theta_vs_Mu.txt')))
    # print(s.join((datapath, '\\', line, '_Fortran_Theta_vs_Mu.txt')))
    # print(s.join((directory, '\\', line, '_Fortran_Theta_vs_Mu.txt')))
    assert len(testdata) == len(outdata)
    assert all(abs(x-y) < 1e-6 for l0, l1 in zip(testdata, outdata)
                               for x, y in zip(l0, l1) )
    exit(process.returncode)
  finally:
    pass

def fixture_test(inpath, execpath, testpath, execrunner=None):
  """ Runs simulation and gets output. """
  from os.path import join, exists
  from os import makedirs, getcwd
  from shutil import rmtree
  from subprocess import Popen, PIPE
  from sys import exit

  directory = join(getcwd(),'zacros_testoutput',testpath)
  exectorun = get_exectorun(execpath, execrunner)
  if not exists(directory): makedirs(directory)

  try: 
    datapath = join(inpath, join('data', testpath))
    copyfiles(datapath, directory)

    process = Popen(exectorun, cwd=directory, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate(input=None)

    print("Stdout from this test was:")
    print(stdout)
    print("Stderr from this test was:")
    print(stderr)
    exit(process.returncode)

  finally: pass

def get_number_of_threads():
  import multiprocessing
  from os import environ
  if (environ.get('OMP_NUM_THREADS')):
    OMP_NUM_THREADS = str(environ.get('OMP_NUM_THREADS'))
  else:
    OMP_NUM_THREADS = str(multiprocessing.cpu_count())
  return OMP_NUM_THREADS

def get_exectorun(execpath, execrunner):
  OMP_NUM_THREADS = get_number_of_threads()
  if execrunner:
    exectorun = [execrunner, '-n', '1', '-d', OMP_NUM_THREADS, execpath]
  else:
    exectorun = [execpath]
  return exectorun


if __name__ == '__main__':
  from sys import argv
  from os.path import dirname, abspath
  import argparse
  parser = argparse.ArgumentParser(description='Run MesoApprox tests.')
  parser.add_argument("execpath", help="Path to MesoApprox executable")
  parser.add_argument("testnumber", help="Test number")
  parser.add_argument("--execrunner",
                      help="MPI runner to execute \
                      tests ie. aprun on Cray systems")
  group = parser.add_mutually_exclusive_group()
  group.add_argument("--fixture", help="Enable fixture testing",
                     action="store_true")
  group.add_argument("--Nevents", type=int, help="Total number of KMC events to\
                      complete in a resume test. Implies resume testing")
  args = parser.parse_args()

  testbasedir = abspath(dirname(argv[0]))
  # There is no good way to get the absolute path of the main script
  # from argparse. So just use this.

  if (args.execrunner):
    execrunner = args.execrunner
  else:
    execrunner = None

  execpath = abspath(args.execpath)
  if (args.Nevents):
    test(testbasedir, execpath, args.testnumber,
         resumetest=True, Neventstarget=args.Nevents, execrunner=execrunner)
  elif(args.fixture):
    fixture_test(testbasedir, execpath, args.testnumber, execrunner=execrunner)
  else:
    test(testbasedir, execpath, args.testnumber, execrunner=execrunner)
