# kuramoto
An efficient C implementation of the standard Kuramoto coupled oscillators system

<img src="formula.png">

with C and Matlab interfaces. Simulations using the Euler and classic Runge-Kutta ("RK4") methods are available. There are C and Matlab programs to demonstrate usage of the interfaces to the **kuramoto** library.

### Building and installation
You will need the [Make](https://www.gnu.org/software/make/) build tool installed on your system, and the Matlab `mex` and `makemex` executables (in the `bin` directory of your Matlab installation) on your system executable path. [^1]

To build and install the kuramoto shared library, in a terminal navigate to the **kuramoto** root directory and run
```
make -C C/lib && sudo make -C C/lib install
```
This builds, and then installs the library. The default installation directory is /usr/local/lib (you will need admin privileges for this).  Alternatively, you may specfify an installation path by setting the PREFIX environmental variable before running `make install`.

To build the C demo, run
```
make -C C/demo
```
and to build the Matlab interface, run
```
make -C matlab
```
### Testing
To test the C interface, try running
```
./kuramoto demo
```
and
```
./kuramoto audio
```
in the **kuramoto**/C/demo directory. To test the Matlab interface, run
```
kuramoto_demo
```
at the Matlab prompt in the **kuramoto**/matlab directory.

### Contributing
If someone would like to contribute a Python (NumPy) interface, please contact the maintainer. This should be reasonably starightforward using  [ctypes](https://docs.python.org/3/library/ctypes.html).

Lionel Barnett: lionelb@sussex.ac.uk
[^1]: Building/installing targets a POSIX environment. Tested with gcc (Linux 64-bit), mingw-w64/MSYS2 (Windows 64-bit) and Clang/LLVM (macOS). Note that MSVC (still) does not fully support the C99 standard, so may be problematic.
