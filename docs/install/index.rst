.. _install:

============
Installation
============

------------------
Build Requirements
------------------

Linux
-----
To build the C++ library on linux, a C++ compiler is required, which supports
the C++17 standard. Any remotely recent version of g++ or clang should fit
the bill for this. You will also need cmake >= 3.9.

In addition, if you would like to build the Python API, you will also need to
ensure that the Python development headers are installed on your system.

Debian / Ubuntu / LinuxMint
+++++++++++++++++++++++++++
To ensure you have all the build requirements on a Debian based system, run the
following command:

.. code-block:: sh

  sudo apt install g++ cmake python3-dev

Fedora / CentOS / RedHat
++++++++++++++++++++++++
To ensure you have all the build requirements on a RedHat based system, run the
following command:

.. code-block:: sh

  sudo dnf install g++ cmake python3-devel

Arch / Manjaro
++++++++++++++++++++++++
To ensure you have all the build requirements on an Arch based system, run the
following command:

.. code-block:: sh

  sudo pacman -S g++ cmake python3

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
  mkdir build && cd build
  cmake -DCMAKE_BUILD_TYPE=Release ..
  make -j
  sudo make install

Build Options
-------------
PapillonNDL has three main build options, which may be activated through the
cmake command. They are listed here:

PNDL_PYTHON
  This is used to build the Python bindings. This is on by default!

PNDL_TESTS
  This is used to build the unit tests, and is turned off by default.

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

-------
Windows
-------
Building PapillonNDL natively on Windows is possible, but not recommended. If
possible, use the Windows subsystem for linux (WSL), and then follow the linux
build instruction.

If for some reason you MUST build it natively for Windows, then you will need
to have Visual Studio 2019 installed. You must also have the Python developement
kit installed, which can be done inside of Visual Studio. When you open the source
directory inside of Visual Studio, it should automatically recognize the cmake
file, and allow you to build the library. I do not use Windows, or VS, but I was
able to build the library this way without a problem inside a Windows VM.
