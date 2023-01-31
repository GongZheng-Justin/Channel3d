# Channel3d
&emsp;**Channel3d** is an efficient second- /fourth-order finite-difference direct numerical simulation (DNS) solver with versatile viscous treatments, also with the ability to handle different boundary conditions:

* Versatile viscous treatments. The viscous term can be handled full implicitly, full explicitly, or partial implicitly
* Have the ability to handle periodic, no-slip and free-slip boundary conditions
* Second-order spatial accuracy in non-periodic directions, and fourth-order scheme is also available for periodic directions 
* FFT-based ([FFTW](https://github.com/FFTW/fftw3) used here) method is used for Pressure Poisson Equation (PPE)
* Alternating Direction Implicit (ADI) method is adopted for Helmholtz equations
* MPI parallelization by means of pencil distributed decomposition, using [2DECOMP&FFT](https://github.com/numericalalgorithmsgroup/2decomp_fft)

## How to cite **Channel3d**?
&emsp;Z. Gong and X. Fu, A pencil distributed direct numerical simulation solver with versatile treatments for viscous term, Computers and Mathematics with Applications, 2021;100:141-151, https://doi.org/10.1016/j.camwa.2021.09.003.

## Notice :newspaper: 
&emsp;**Channel3d** has been integrated into the CFD-DEM sovler [**CP3d**](https://github.com/GongZheng-Justin/CP3d) as a sub-solver.

## Installation :briefcase:
&emsp;During developing this solver, I often try my best to make it easy-to-understand and easy-to-use. As for compilation, present solver only has the following two prerequisites:

* MPI
* Gfortran/Intel Fortran (Supporting Fortran 2003 or higher version)

&emsp;**FFTW-3.3.9** library has been explicitly included in the directory `./src/ThirdParty/fftw/`, but it is strongly recommended to recomplied FFTW-3.3.9 for the first use. After entering the folder `Channel3d-master/` in terminal, you can compile the code as follows:
```
1. chmod a+x ./mymake.sh
2. ./mymake.sh
3. choose the correct compiler you use, and the executable you want to compile, following guidances printed in the terminal
```
&emsp;Yon can also compile the `interpolateField` code in the folder `./Tool/interpolateField/` by typing:
```
1. cd ./Tool/interpolateField
2. chmod a+x ./makeInterp.sh
3. ./makeInterp.sh
4. choose the correct compiler you use, and the executable you want to compile, following guidances printed in the terminal
5. cd ../..
```
&emsp;If the compiling process successfully, the executable file(s) `channel2nd`/`channel4th` will be appeared in the folder `Channe3d-master/`, and `interpolateField` will be included in the folder `./Tool/interpolateField/`.

## Usage :book:
&emsp;After compiling the code successfully, you can run the executable file like that:
```
mpirun -n [np] [exeName] [inputFile]
```
&emsp;Here:
* `np` denotes the number of processors you use
* `exeName` stands for specific executable file name, namely `channel2nd` or `channel4th`
* `inputFile` is the name string for the input parameter file  

&emsp;For instance, if you want to run the lid-driven cavity case, you can type the following words in your terminal:
```
mpirun -n 4 ./channel2nd ./Input/LidDrivenCavity.prm
```
### Input file
&emsp;The input file examples are stored in the folder `./Input/`. See `./doc/channel2nd_prm.md` and `./doc/channel4th_prm.md` for detailed descriptions to the input file for second-order and fourth-order scheme respectively.

### A complete example
&emsp;See `./doc/Toturial_for_channel_turbulence_4th.pdf` for a complete example to the wall-bounded turbulence at Re_tau =180 using fourth-order scheme, including the postprocessing, and verification.

## To do list :muscle:

* Fourth-order scheme options for non-periodic directions
* Hybrid MPI/OpenMP parallelization and GPU acceleration  
* **Channel3d** presented here is a part of my integral project: [**CP3d**](https://github.com/GongZheng-Justin/CP3d), Channnel-Particle 3d, which is still under development, and will be open-source in the coming future. My ultimate aim for **CP3d** project is to develop an efficient and easy-to-use  particle-laden flow solver, including one-way, two-way, and full four-way coupling methods, by combination with Basset history force model, discrete element method (DEM), and immersed boundary method (IBM).

## Contact and Feedback :email:
&emsp;If you have any question, or want to contribute to the code, please don't hesitate to contact me: Zheng Gong (gongzheng_justin@outlook.com)

## Showcase
&emsp;Following picture shows the transition process from laminar to turbulent flow for open-channel case (Re_tau =180). The top half of the picture presents the normalized streamwise velocity field, and the bottom half shows the isosurface for Q vortex criterion value at Q+ = 1/80, colored by streamwise velocity.
![](doc/cha180_small.gif)
