# Notes from Hunter
I found these scripts on the IAEA website [here](https://www-nds.iaea.org/TSL_LibGen/). They
are written by Jose Ignacio Marquez Damian, who produced the new evaluations for H2O and
D2O in ENDF/B-VIII. They are used to perform interpolation of the phonon spectra and
oscillators in order to produce thermal scattering laws for H in H2O and D D2O at
temperatures for which a phonon spectrum isn't provided. It was originally written in Python2,
but I have converted it to run with Python3, and it appears to be working just fine. The
original script `leapr_interpolator.py` could only produce a LEAPR input file for a single
temperature. I have made a modified version, `new_leapr_interpolator.py`, which is able
to produce a LEAPR input file for multiple temperatures at a time.

# Original Readme Text

Files contained in this directory:

- Directories:
* original_inputs/
  Original input files
 
* xml/
  Processed inputs in XML format
 
- Utility programs:
* leapr2xml.py
 Converts LEAPR inputs into XML files for internal usage.
Example:
  ./leapr2xml.py original_inputs/hh2o-endf6.njoy test.xml

* leapr_check.py
 Checks a LEAPR input in two steps: first it converts it to XML, then it checks for logical errors.
Example:
  ./leapr_check.py original_inputs/hh2o-endf6.njoy
 
* leapr_interpolator.py
 Interpolates a LEAPR model in XML format for a given temperature and produces a LEAPR input file.
 The program checks the original model and the interpolated model.
 ./leapr_interpolator.py xml/hh2o-endf6.xml test.leapr 301
 
* run_leapr.py
 Creates an interpolated LEAPR input file and runs it locally or in a remote server.
 
 The following variables need to be modified in the Python file:
 
 server: points to the server to run NJOY. A blank string means NJOY will run locally.
 tmpdir: points to the temporary directory in the server. "/tmp" is probably good.
 njoy_exec: points to the NJOY executable in the server.
 dryrun: if set to True, creates the input file but does not run it. Set it to False to actually run NJOY.

