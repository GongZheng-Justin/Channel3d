# Channel3d
&emsp;**Channel3d** is an efficient second- /fourth-order finite-difference direct numerical simulation (DNS) solver with versatile viscous treatments, also with the ability to handle different boundary conditions:

* Versatile viscous treatmens. The viscous term can be handled full implicitly, full explicitly, or partial implicitly
* Have the ability to handle periodic, no-slip and free-slip boundary conditions
* Second-order spatial accuray in non-periodic directions, fourth-order scheme is also available for periodic directions 
* FFT-based ([FFTW](https://github.com/FFTW/fftw3) used here) method is used for Pressure Poisson Equation (PPE)
* Alternating direction implicit (ADI) is adopted for Helmholtz equations
* MPI parallelization by means of pencil distributed decomposition, using [2DECOMP&FFT](http://www.2decomp.org/)

## Installation :briefcase:
&emsp;During developing this solver, I often try my best to make it easy-to-understand and easy-to-use. As for compilation, present solver only has the following two prerequisites:

* MPI
* Gfortran/Intel Fortran (Supporting Fortran 2003 or higher version)

&emsp;**FFTW-3.3.8** library has been explicitly included in the direction `./src/ThirdParty/fftw/`, so compiling and additional linking to external FFTW are avoided. After entering the folder `Channe3d-master` in terminal, you can compile the code as follows:
```
1. chmod a+x ./mymake.sh
2. ./mymake.sh
3. choose the correct compiler you use, and the executable you want to compile, following guidances printed in the terminal
```
## Usage :book:
&emsp;

## To do list :muscle:

* Fourth-order scheme options for non-periodic directions
* Adding a passive scalar transport solver
* Hybrid MPI/OpenMP parallelization and GPU acceleration  
* **Channel3d** presented here is a part of my integral project: **CP3d**, Channnel-Particle 3d, which is still under development, and will be open source in the coming future. My ultimate aim for **CP3d** project is to develop an efficient and easy-to-use  particle-laden flow solver, including one-way, two-way, and full four-way coupling methods, by combination with Basset history force model, discrete element method (DEM), and immersed boundary method (IBM).

## Acknowledgements :clap:
&emsp;Since Sep 2019, when I finally decided to develop my own CFD-DEM code from scratch, I have learnt quite a lot from the following really kind researchers (**in alphabetical sequence**):

* [Dr. Costa](https://p-costa.github.io/) from University of Iceland, and his second-order DNS code [CaNS](https://github.com/p-costa/CaNS), also his papers on IBM approach.
* [Dr. He](https://www.engineering.iastate.edu/people/profile/phe/) from Iowa State University, and his fourth-order DNS solver [HercuLES](https://github.com/friedenhe/hercules).
* [Prof. Ji](http://faculty.tju.edu.cn/ChunningJi/en/index.htm) from Tianji University, on the fruitful discussion about the particle IBM method, and on the accsee to their in-house DNS/LES-Solid interaction code **_cgLES_**.
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
