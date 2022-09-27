# kuramoto
An efficient C implementation of the standard Kuramoto coupled oscillators system

<img src="formula.png">

with Matlab interface. Simulations using the Euler and classic Runge-Kutta ("RK4") methods are available.

### Building
You will need the [Make](https://www.gnu.org/software/make/) build tool installed on your system, and the Matlab `mex` and `makemex` executables (in the 'bin' directory of your Matlab installation) on your system executable path. Then in a terminal, navigate to the **kuramoto** installation directory and run `make`.[^1]

### Installation
Building creates a shared library (libkuramoto.so in Linux and macOS, libkuramoto.dll in Windows) in the 'lib' sub-directory. You will need to ensure that this directory is on your dynamic link path: this is specified by the environmental variable LD_LIBRARY_PATH in Linux, DYLD_LIBRARY_PATH in macOS and just PATH in Windows (consult your OS documentation about setting environmental variables). Alternatively, you can copy the shared library to a standard system library directory (you may need admin privileges for this)[^2].

If you do not wish to bother with installation, a quick fix is to run Matlab, e.g., on Linux as

> LD_LIBRARY_PATH=<kuramoto directory>/lib:$LD_LIBRARY_PATH matlab

### Running
To test the Matlab interface, we recommend that you run the `kuramoto_demo.m` script in the kuramoto/Matlab directory.

### Contributing
If someone would like to contribute a Python (NumPy) interface, please contact the maintainer.

Lionel Barnett: lionelb@sussex.ac.uk
[^1]: Tested with gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and Clang/LLVM (macOS). Note that MSVC (still) does not fully support the C99 standard, so may be problematic.
[^2]: On Windows you may have to register the dll. I am not familiar with Windows installation requirements.
