# Papillon Nuclear Data Library
[![Documentation Status](https://readthedocs.org/projects/papillon-ndl/badge/?version=latest)](https://papillon-ndl.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/License-CeCILL%20v2.1-brightgreen)](https://cecill.info/index.en.html)

The Papillon Nuclear Data Library (NDL) is a C++/Python library for reading,
sampling, and interacting with continuous energy nuclear data, stored in the ACE
nuclear data format.

For examples of how to use the library in both C++ and Python, take a look at
the [documentation site](https://papillon-ndl.readthedocs.io). That is where you
will also be able to fined more detailed installation instructions.

## License
PapillonNDL is provided under the terms and conditions of the CeCILLv2.1
license. This is a French equivalent of the GPLv3 license, and is explicitly
compatible with both the GPLv2 and v3. The French version of this license is
provided in the LICENSE file, along with the equally valid English version, which
may be found in the LICENSE-ENGLISH file. The CeCILLv2.1 is approved by the FSF.
More information about this license may be found [here](https://cecill.info/).

## Dependencies
To build and install the library a Unix-like operating system with cmake >= 3.11
is required, along with a C++ compiler which supports the C++17 standard. The
recommended compiler is Clang >= 6 or GCC >= 7 should suffice. In order to
build the Python interface, Python >= 3.5 should be installed on your system, in
addition to the Python development libraries and header files.

There are two posibilities for building the library on Windows. The first is to use
the Windows Subsytem for Linux (WSL), where the library may be installed by following
the simple Linux build instructions, typically without problem. The second option
is to build the library to run natively on Windows. This requires having Visual
Studio installed to compile the library. In addition, if you want to build the
Python bindings for Windows, you need to ensure the Python development kit has been
installed in Visual Studio as well.

Tests are not built by default, and should only be needed for developers. You
can turn them on by using ```-DPNDL_TESTS=ON``` with cmake.

## Install
To build PapillonNDL on a Unix-like system, navigate to the directory where you
would like to keep the source files, and then run the following commands:
```
$ git clone https://github.com/HunterBelanger/papillon-ndl.git
$ cd papillon-ndl
$ mkdir build && cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ sudo make install
```
This will install the libraries and header files to the UNIX default locations
in ```/usr/local/```.

If you do NOT want to build the Python bindings for PapillonNDL, then you should
add the flag ```-DPNDL_PYTHON=OFF``` to the cmake command.
