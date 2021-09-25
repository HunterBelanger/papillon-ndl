.. _install:

============
Installation
============

------------------
Build Requirements
------------------

To build the C++ library on linux, a C++ compiler is required, which supports
the C++20 standard. On Linux systems, the only current option which satisfied
this is GCC >= 11. On Windows, you will need MSVC >= 19.29. You will also need
cmake >= 3.11, whether you use Windows or linux.

In addition, if you would like to build the Python API, you will also need to
ensure that the Python development headers are installed on your system.

Fedora 35+
------------------------
To ensure you have all the build requirements on a RedHat based system, run the
following command:

.. code-block:: sh

  sudo dnf install g++ cmake python3-devel

Arch / Manjaro 21+
------------------------
To ensure you have all the build requirements on an Arch based system, run the
following command:

.. code-block:: sh

  sudo pacman -S g++ cmake python3

Other Linux Distributions
-------------------------
Unfortunately, GCC-11 is quite new, and many distrubtions do not yet ship this
compiler toolchain, or have specific installation commands. If you are not
using one of hte distributions listed above, first check to see if there is
a specific method to install GCC-11 form your distributions repositories. If
not, I recommend building GCC-11 from source. This is outside the scope of these
installation instruction, but instructions are readily found online.

Windows
--------
There are no quick and easy commands to install the necessary dependencies on
Windows sadly. You should insall Visual Studio 2019, and when performing the
installation, be sure to also install the Python development libraries. You
will also have to download and install cmake separately as well. If you already
have Visual Studio 2019 installed, you should make sure that it has been
updated to the most recent version, or that the compiler version is at least
MSVC >= 19.29.

------------------
Getting the Source
------------------
The source can either be downloaded directly from GitHub, or via git. If using
git, first navigate to the directory where you would like to keep the source
folder, and then run the command:

.. code-block:: sh

  git clone https://github.com/HunterBelanger/papillon-ndl.git

This will download the soruce for the latest development version. If you would
prefer to use the latest stable version, then subsequently run

.. code-block:: sh

  cd papillon-ndl
  git checkout master

If you would like to install a particular version of the library, you can also
checkout the desired version tag:

.. code-block:: sh

  git checkout v0.1.1

-----
Build
-----

Quick'n Dirty
-------------
If you are looking to build the C++ and Python libraries, and don't want to do
any tinkering, then just run these commands, and you should be good to go !

.. code-block:: sh

  cd papillon-ndl
  cmake -E make_directory build
  cd build
  cmake -DCMAKE_BUILD_TYPE=Release ..
  cmake --build . --target install

On linux systems, you may need to add ``sudo`` to the beginning of the command,
to allow cmake to install the files into the ``/usr/local`` directory.

Build Options
-------------
PapillonNDL has three main build options, which may be activated through the
cmake command. They are listed here:

PNDL_PYTHON
  This is used to build the Python bindings. This is on by default!

PNDL_TESTS
  This is used to build the unit tests, and is turned off by default.

PNDL_SHARED
  Builds a shared library, as opposed to a static library. This is turned on by
  default. When building on Windows, this will automatically be turned off.

PNDL_INSTALL
  Adds and exports the installation targets for PapillonNDL. This is truned on
  by default, but can be turned off if using PapillonNDL as a build dependency
  in another project.

Several other standard cmake options will also be usefull in many cases, and
are therefore listed here:

CMAKE_BUILD_TYPE
  If you are looking for sane optimizations (``-O2``) in a normal build, set
  this to ``Release``. When doing development, it is often adventageous to set
  this to ``Debug``, which provides debug symbols.

CMAKE_INSTALL_PREFIX
  This is the location to where the libararies and header files will be
  installed on your system. On linux systems, this is usually ``/usr/local``.
  If you want to change it, you can set it with this command.

As an example, if we wanted to build PapillonNDL without the Python bindings,
in debug mode, and install it to our home directory, then when running cmake
we should use:

.. code-block:: sh

  cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/home/user -DPNDL_PYTHON=Off ..

In this case, the header files will be installed to the directory 
``/home/user/include/PapillonNDL``, the C++ library will be installed to
``/home/user/lib``, and the Python library will be installed to a directory similar
to ``/home/user/lib64/python3.9/site-packages`` (the Python version might be
different however).

.. warning::
  
  Make sure the the directory where the Python library was installed is in your
  ``PYTHONPATH`` environement variable ! If it isn't, Python won't be able to
  find the library ! If you open python in your terminal, and can run
  ``import pyPapillonNDL`` without an error, you should be fine. If you do get
  an error, add the path to the ``pyPapillonNDL`` library to your ``PYTHONPATH``.

  This should only be a problem if you used the ``CMAKE_INSTALL_PREFIX`` option to
  install to a different location than the default.
