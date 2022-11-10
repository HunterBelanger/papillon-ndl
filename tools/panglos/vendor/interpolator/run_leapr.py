#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import sys
import string
import random
import leapr_interpolator

tfreeze = 273.15 # K
tcritical = 647 # K

#server = "localhost"
server = "" # if server is empty, runs locally. If not, it runs by ssh.
njoy_exec="./njoy/njoy"
tmpdir = "/tmp/"
dbg = 2 # 0=quiet, 1=brief, 2=detailed
dryrun = False # if True, does not execute NJOY

available_models = {"HH2O-ENDF6":"xml/hh2o-endf6.xml", "HH2O-ENDF7":"xml/hh2o-endf7.xml", "HH2O-ENDF8":"xml/hh2o-cab.xml", "DD2O-ENDF8":"xml/dd2o-cab.xml", "OD2O-ENDF8":"xml/od2o-cab.xml", "DD2O-ENDF7":"xml/dd2o-endf7.xml"}

def random_string_generator(size=8, chars=string.ascii_uppercase + string.digits):
  return ''.join(random.choice(chars) for x in range(size))

def run_input(basename, nout):
#
# Create directory with random name in running server
#
    rundir = tmpdir+random_string_generator()

    leapr_input = basename + ".leapr"
    njoy_input = basename + ".njoy"
    out_file = basename + ".out"
    status_file = basename + ".status"
    endf_file = basename + ".endf"
    
    f = open(njoy_input, "w")
    f.write("#!/bin/env bash\n")
    f.write("#\n")
    f.write("# Bash script to create the input and run njoy\n")
    f.write("#\n")
    f.write("NJOY=%s\n" % (njoy_exec))
    f.write("cat>input<<EOF\n")
    f.write(open(leapr_input, "r").read())
    f.write("EOF\n")
    f.write("${NJOY} < input\n")
    f.write("cp tape%s %s\n" % (nout, endf_file))
    f.write("if [ $? -eq 0 ]\n")
    f.write("then\n")
    f.write("  echo -n 'OK' > %s\n" % status_file)
    f.write("else\n")
    f.write("  echo -n 'NOT OK' > %s\n" % status_file)
    f.write("fi\n")
    f.close()
    
    if (server == ""):
      command = "mkdir %s > /dev/null 2>&1" % (rundir)
      if (dbg>=1):
        print ">>> Creating random directory %s" %(rundir)
    else:
      command = "ssh -q %s 'mkdir %s' > /dev/null 2>&1" % (server, rundir)
      if (dbg>=1):
        print ">>> Creating random directory %s in %s" %(rundir, server)
    if (dbg>=2):
      print command
    if not dryrun:
      os.system(command)
#
# Copy input file
#
    if (server == ""):
      command = "cp  %s %s/ > /dev/null 2>&1" % (njoy_input, rundir)
    else:
      command = "scp -q %s %s:%s/ > /dev/null 2>&1" % (njoy_input, server, rundir)
    if (dbg>=1):
      print ">>> Copying input file"
    if (dbg>=2):
      print command
    if not dryrun:
      os.system(command)
#
# Run input file
#
    if (server == ""):
      command = "cd %s; sh ./%s > %s 2>&1" % (rundir, njoy_input, out_file)
      if (dbg>=1):
        print ">>> Running NJOY"
    else:
      command = "ssh -q %s 'cd %s; sh ./%s > %s 2>&1'" % (server, rundir, njoy_input, out_file)
      if (dbg>=1):
        print ">>> Running NJOY in %s" % server
    if (dbg>=2):
      print command
    if not dryrun:
      os.system(command)
#
# Copy ENDF-6 file back
#
    if (server == ""):
      command1 = "cp %s/%s . > /dev/null 2>&1" % (rundir, endf_file)
      command2 = "cp %s/%s . > /dev/null 2>&1" % (rundir, out_file)
      command3 = "cp %s/%s . > /dev/null 2>&1" % (rundir, status_file)
    else:
      command1 = "scp -q %s:%s/%s . > /dev/null 2>&1" % (server, rundir, endf_file)
      command2 = "scp -q %s:%s/%s . > /dev/null 2>&1" % (server, rundir, out_file)
      command3 = "scp -q %s:%s/%s . > /dev/null 2>&1" % (server, rundir, status_file)
    if (dbg>0):
      print ">>> Retrieving files"
    if (dbg>=2):
      print command1
      print command2
      print command3
    if not dryrun:
      os.system(command1)
      os.system(command2)
      os.system(command3)
#
# Check run
#
    if not dryrun:
      status = open(status_file, "r").read()
    else:
      status = "OK"
#
# Delete run directory and temporary files
#
    if (server == ""):
      command1 = "rm -rf %s > /dev/null 2>&1" % (rundir)
      if (dbg>=1):
        print ">>> Erasing temporary directory %s/" % (rundir)
    else:
      command1 = "ssh -q %s 'rm -rf %s' > /dev/null 2>&1" % (server, rundir)
      if (dbg>=1):
        print ">>> Erasing temporary directory %s/ in %s" % (rundir, server)
    command2 = "rm  %s %s > /dev/null 2>&1" % (status_file, leapr_input)
    if (dbg>=2):
      print command1
      print command2
    if not dryrun:
      os.system(command1)
    os.system(command2)

    return status
    
if __name__ == '__main__':
  if "server" not in globals():
    print "Configuration error"
    print "Edit %s and define the variable \"server\"" %(sys.argv[0])
    raise SystemExit
  if "njoy_exec" not in globals():
    print "Configuration error"
    print "Edit %s and define the variable \"njoy_exec\" to point to the NJOY executable in the server" %(sys.argv[0])
    raise SystemExit
  if "tmpdir" not in globals():
    print "Configuration error"
    print "Edit %s and define the variable \"tmpdir\" to point to the temporary directory in the server" %(sys.argv[0])
    raise SystemExit
  if (len(sys.argv)<=2):
    print "Creates a LEAPR input file at a given temperature and runs it in a remote server"
    print "Usage: %s model temperature" %(sys.argv[0])
    print "Model has to be one of: %s" % (", ".join(available_models.keys()))
    raise SystemExit
  model_choice = sys.argv[1]
  model_choice = model_choice.upper()
  if (not (model_choice in available_models.keys())):
    print "Creates a LEAPR input file at a given temperature and runs it in a remote server"
    print "Usage: %s model temperature" %(sys.argv[0])
    print "Model has to be one of: %s" % (", ".join(available_models.keys()))
    raise SystemExit
  try:
    temp = float(sys.argv[2]) # K
    temp = round(temp,2)
  except ValueError:
    print "Error: ", sys.argv[2], " is not a real number"
    raise SystemExit

  if dryrun:
    print ""
    print "WARNING: this is just a test run. To run files in the server, edit %s and set the variable \"dryrun\" to False" % (sys.argv[0])
    print ""
  
  if ((temp < tfreeze) or (temp > tcritical)):
    print "Error: temperature outside of range for liquid water"
    raise SystemExit
  
  cur_dir = os.path.dirname(os.path.realpath(__file__))
  if server == "":
    njoy_exec = os.path.join(cur_dir, njoy_exec)
    if not os.path.isfile(njoy_exec):
      print "Error: NJOY executable not found at %s" % njoy_exec
      print "Follow the instructions at https://njoy.github.io/Build/index.html to obtain NJOY2016"
      raise SystemExit

  basename = "%s-%i" % (model_choice,int(round(temp*100,0)))
  leapr_input = basename + ".leapr"
  model = available_models[model_choice]
  model = os.path.join(cur_dir, model)
  nout = leapr_interpolator.getnout(model)

  if (dbg >= 1):
    print ">>> Creating the input by interpolation"
  leapr_interpolator.leapr_interpolator(model, leapr_input, temp)

  status = run_input(basename, nout)

  if (dbg >= 1):
    print ">>> The status of the run is: %s" % status

