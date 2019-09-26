# Linear Constrained Minimization Solver

This is an interior point method for linearly constrained minimization.

## Pre-requisites

1. You need to create a `./lib` folder and include the following static libraries
   * `libalgencan.a` generated from [Algencan](https://www.ime.usp.br/~egbirgin/tango/downloads.php)
   * `libhsl_ma48.a` generated form [HSL MA48](http://www.hsl.rl.ac.uk/catalogue/hsl_ma48.html)
   * `libhsl_ma57.a` generated from [HSL MA57](http://www.hsl.rl.ac.uk/catalogue/hsl_ma57.html)
   * `libmetis.a` generated from [MeTIS library](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz)
   
2. You need to create a `./include` folder and include
   * `hsl_ma48_double.mod` (generated after compiling HS MA48)
   * `hsl_ma57_double.mod` (generated after compiling HS MA57)
   * `hsl_zd11_double.mod` (generated after compiling HS MA48)
   
3. You need to include [mc58ad.f](http://www.hsl.rl.ac.uk/catalogue/mc58.html) under `./sources/hsl`
