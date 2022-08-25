# kuramoto
An efficient C implementation of the standard Kuramoto coupled oscillators system

<img src="formula.png">

with Matlab interface. Simulations using the Euler and classic Runge-Kutta ("RK4") methods are available.

### Building
You will need the [Make](https://www.gnu.org/software/make/) build tool installed on your system, and the Matlab `mex` and `makemex` executables (in the 'bin' directory of your Matlab installation) on your system executable path. Then in a terminal, navigate to the **kuramoto** installation directory and run `make`.[^1]

### Running
To test the Matlab interface, we recommend that you run the `kuramoto_demo.m` script in the kuramoto/Matlab directory.

### Contributing
If someone would like to contribute a Python (NumPy) interface, please contact the maintainer.

Lionel Barnett: lionelb@sussex.ac.uk
[^1]: Tested with gcc (Linux 64-bit), mingw-w64 (Windows 64-bit) and Clang/LLVM (macOS). Note that MSVC (still) does not fully support the C99 standard, so may be problematic.
