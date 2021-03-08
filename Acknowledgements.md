Since Oct 2019, when I finally decided to develop my own CFD-DEM code from scratch, I have learnt quite a lot from the following really kind researchers and helpful references:
1. [Prof. Paolo Orlandi](http://dma.ing.uniroma1.it/users/orlandi/resume.html) from Sapienza University of Rome, and his book ***Fluid flow phenomena: a numerical toolkit***, besides the [attached CFD code](http://dma.ing.uniroma1.it/users/orlandi/diskette.tar.gz).
2. [Dr. hamid reza norouzi](https://www.researchgate.net/profile/Hamid-Norourzi) from University of Tehran, and his book ***Coupled CFD‐DEM Modeling: Formulation, Implementation and Applimation to Multiphase Flows***, besides the [attached DEM code](https://www.wiley.com//legacy/wileychi/norouzi/form.html?type=SupplementaryMaterial).
3. [Incompact3d](https://github.com/xcompact3d/Incompact3d)
4. [AFid](https://github.com/PhysicsofFluids/AFiD)
5. [CaNS](https://github.com/p-costa/CaNS)
6. [HercuLES]()


ThirdParty part:

  1. The 'decomp2d' in directory './ThirdParty/' was copied from the open source code 'Incompact3d-3.0'
            ( Downloaded from: https://github.com/xcompact3d/Incompact3d, about Feb 2020 )

     Some necessary modifications were done:
        (1) The 'fft-*' parts are deleted.

        (2) The files 'mem_merge.f90' and 'mem_split.f90' are deleted.

        (3) Add two files:  'transpose_x_to_z.inc', and 'transpose_z_to_x.inc' from the open source code 'cans'
            ( Downloaded from: https://github.com/p-costa/CaNS,  April 17, 2020 )
            And this two files are almost the same with another open source code 'AFiD', except for the data type.
            ( Downloaded from: https://github.com/PhysicsofFluids/AFiD, April 26, 2020 )

        (4) Some new lines are added into the 'decomp_2d.f90' to make the codes can handle 'transpose_x_to_z' and 'transpose_z_to_x' directly( According to 'cans' and 'AFiD') . 

        (5) I added a new 'mydecomp_2d_extra.f90' to handle the update_halo more conveniently.


  2. The 'fftw3' in directory './ThirdParty/' is a version of fftw-3.3.8. It was added into the 'src' to avoid linking to external 'include' and 'lib'.        
   

    Gong Zheng
    2021/03/08
