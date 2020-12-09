# Papillon Nuclear Data Library
![License](https://img.shields.io/badge/License-CeCILL%20v2.1-brightgreen)

The Papillon Nuclear Data Library (NDL) is used for reading, sampling,
and interacting with continuous energy nuclear data, stored in the ACE
nuclear data format.

### Sampling
To sample the angle and energy distributions, an instance of an
```std::function<double()>``` must be passed to the distributions.
This function pointer is the library's access to a pseudo random
number generator. Therefore, if ```rng``` is such an instance, ```rng``` must
meet the requirments that any call of ```rng()``` returns a random double over
the interval [0,1), and it must be possible to call ```rng()``` an indefinite
number of times. If the passed object does not meet these requirments, the
behaviour of the library is undefined. Using this mechanism allows greater
flexibility to the user, in choosing the random number generator of their
choice, and allows the library to be completely disconnected from the random
number generation process, allowing the distributions to know how to sample
themselves, but not know or assume anything about random number generation.

## Dependencies
To build and install the library, cmake >= 3.9 is required, along with a C++
compiler which supports the C++17 standard. The recommended compiler is
Clang >= 5, though GCC >= 6 should suffice. Inorder to build the Python
interface, Python >= 3.5 should be installed on your system as well.

Building the unit tests (using ```-DPNDL_TESTS=ON``` when calling cmake)
requires that Google test already be installed on the system, and is not
provided. Tests are not built by default, and should only be needed for
developers, therefore this is not required for most users.

## Install
To build PapillonNDL, navigate to the directory where you would like to keep
the source files, and then run the following commands:
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

## License
PapillonNDL is provided under the terms and conditions of the CeCILLv2.1
license. This is a French equivalent of the GPLv3 license, and is explicitly
compatible with both the GPLv2 and v3. The French version of this license is
provided in the LICENSE file, along with the equaly valid English version, which
may be found in the LICENSE-ENGLISH file. The CeCILLv2.1 is approved by the FSF.
More information about this license may be found [here](https://cecill.info/).
