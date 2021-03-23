# Channel3d
&emsp;**Channel3d** is an efficient second- /fourth-order finite-difference direct numerical simulation (DNS) solver with versatile viscous treatments, also with the ability to handle different boundary conditions:

* Versatile viscous treatments. The viscous term can be handled full implicitly, full explicitly, or partial implicitly
* Have the ability to handle periodic, no-slip and free-slip boundary conditions
* Second-order spatial accuray in non-periodic directions, fourth-order scheme is also available for periodic directions 
* FFT-based ([FFTW](https://github.com/FFTW/fftw3) used here) method is used for Pressure Poisson Equation (PPE)
* Alternating Direction Implicit (ADI) method is adopted for Helmholtz equations
* MPI parallelization by means of pencil distributed decomposition, using [2DECOMP&FFT](http://www.2decomp.org/)
![](doc/cha180_small.gif)

## Installation :briefcase:
&emsp;During developing this solver, I often try my best to make it easy-to-understand and easy-to-use. As for compilation, present solver only has the following two prerequisites:

* MPI
* Gfortran/Intel Fortran (Supporting Fortran 2003 or higher version)

&emsp;**FFTW-3.3.8** library has been explicitly included in the directory `./src/ThirdParty/fftw/`, so compiling and additional linking to external FFTW are avoided. After entering the folder `Channe3d-master/` in terminal, you can compile the code as follows:
```
1. chmod a+x ./mymake.sh
2. ./mymake.sh
3. choose the correct compiler you use, and the executable you want to compile, following guidances printed in the terminal
```
&emsp;Yon can also compile the `interploteFiled` code in the folder `./Tool/interploteFiled/` by typing:
```
1. cd ./Tool/interploteFiled
2. chmod a+x ./makeInterp.sh
3. ./makeInterp.sh
4. choose the correct compiler you use, and the executable you want to compile, following guidances printed in the terminal
5. cd ../..
```
&emsp;If the compiling process successfully, the executable file(s) `channel2nd`/`channel4th` will be contained in the folder `Channe3d-master`, and `interpolateField` will be included in the folder `./Tool/interploteFiled/`.

## Usage :book:
&emsp;After compiling the code successfully, you can run the executable file like that:
```
mpirun -n [np] [exeName] [inputFile]
```
&emsp;Here:
* `np` denotes the number of processors you use
* `exeName` stands for specific executable file name, namely `channel2nd` or `channel4th`
* `inputFile` is the name string for the input parameter file  

&emsp;For instatance, if you want to run the lid-driven cavity case, you can type the following words in your terminal:
```
mpirun -n 4 ./channel2nd ./Input/LidCavity.prm
```
### Input file
&emsp;The input file examples are stored in the folder `./Input/`. See `./Input/channel2nd_prm.md` and `./Input/channel4th_prm.md` for detailed descriptions to the input file for second-order and fourth-order scheme respectively.

### A complete example
&emsp;See `./doc/Toturial_for_channel_turbulence_4th` for a complete example to the wall-bounded turbulece at $\sum_{i=0}^N{X_i}$

## To do list :muscle:

* Fourth-order scheme options for non-periodic directions
* Adding a passive scalar transport solver
* Hybrid MPI/OpenMP parallelization and GPU acceleration  
* **Channel3d** presented here is a part of my integral project: **CP3d**, Channnel-Particle 3d, which is still under development, and will be open source in the coming future. My ultimate aim for **CP3d** project is to develop an efficient and easy-to-use  particle-laden flow solver, including one-way, two-way, and full four-way coupling methods, by combination with Basset history force model, discrete element method (DEM), and immersed boundary method (IBM).

## Acknowledgements :clap:
&emsp;Since Sep 2019, when I finally decided to develop my own CFD-DEM code from scratch, I have learnt quite a lot from the following really kind researchers (**in alphabetical sequence**):

* [Dr. Costa](https://p-costa.github.io/) from University of Iceland, and his second-order DNS code [CaNS](https://github.com/p-costa/CaNS), also his papers on IBM approach.
* [Dr. He](https://www.engineering.iastate.edu/people/profile/phe/) from Iowa State University, and his fourth-order DNS solver [HercuLES](https://github.com/friedenhe/hercules).
* [Prof. Ji](http://faculty.tju.edu.cn/ChunningJi/en/index.htm) from Tianji University, on the fruitful discussion about the particle IBM method, and on the access to their in-house DNS/LES-Solid interaction code **_cgLES_**.
* [Dr. Laizet](http://www.imperial.ac.uk/people/s.laizet) from Imperial College London, and their compact FD code [Incompact3d](https://github.com/xcompact3d/Incompact3d).
* [Prof. Marchioli](http://158.110.32.35/) from University of Udine, on the fruitful and continuous discussion about one-way CFD-Particle coupling benchmark and on the access to their [benckmark data](http://158.110.32.35/download/DNS-TEST-CASE/).
* [Dr. Norouzi](https://www.researchgate.net/profile/Hamid-Norourzi) from University of Tehran, and his book **_Coupled CFD‐DEM Modeling: Formulation, Implementation and Applimation to Multiphase Flows_**, besides the [attached DEM code](https://www.wiley.com//legacy/wileychi/norouzi/form.html?type=SupplementaryMaterial).
* [Prof. Orlandi](http://dma.ing.uniroma1.it/users/orlandi/resume.html) from Sapienza University of Rome, and his book **_Fluid flow phenomena: a numerical toolkit_**, besides the [attached CFD code](http://dma.ing.uniroma1.it/users/orlandi/diskette.tar.gz).
* [Dr. Tschisgale](https://www.researchgate.net/profile/Silvio-Tschisgale) from Institute of Air Handling and Refrigeration, on the fruitful and continuous discussion about their IBM approach.
* [Prof. Zhao](http://www.hy.tsinghua.edu.cn/info/1154/1829.htm) from Tsinghua university, on the one-way CFD-Particle coupling benchmark.
* ......

&emsp;Without those researchers' help, I might do nothing but sleep in the dormitory all the days!!!:joy::joy::joy:   
&emsp;Thanks so much again and again !!!

## Contact and Feedback :email:
&emsp;If you have any question, or want to contribute to the code, please don't hesitate to contact me: Zheng Gong (gongzheng_justin@outlook.com)
