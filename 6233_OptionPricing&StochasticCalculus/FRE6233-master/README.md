# FRE6233

A library for financial engineering.

This libray has platform-independent C++ header files (starting with `fms_`) for option pricing and greeks. The corresponding files starting with `xll_` implement Excel add-ins to call the C++ code.

## Build

After cloning the repository and opening the solution, pressing `F5` should build the add-in and start Excel in the debugger with the add-in loaded. You must specify the `x86` platform if you use 32-bit Excel.
