# ThirdParty part
The 'decomp2d' in directory './ThirdParty/' was copied from the open source code [Incompact3d-3.0](https://github.com/xcompact3d/Incompact3d),downloaded in Feb 2020.
Some necessary modifications were done:
* The 'fft-*' parts are deleted.
* The files 'mem_merge.f90' and 'mem_split.f90' are deleted.
* Add two files:  'transpose_x_to_z.inc', and 'transpose_z_to_x.inc' from the open source code [cans](https://github.com/p-costa/CaNS), downloaded in  April 2020. And this two files are almost the same with another open source code [AFiD](https://github.com/PhysicsofFluids/AFiD), downlaoded in April 2020, except for the data type.
* Some new lines are added into the 'decomp_2d.f90' to make the codes can handle 'transpose_x_to_z' and 'transpose_z_to_x' directly( According to 'cans' and 'AFiD') . 
* I added a new 'mydecomp_2d_extra.f90' to handle the update_halo more conveniently.

The 'fftw3' in directory './ThirdParty/' is a version of fftw-3.3.8. It was added into the 'src' to avoid linking to external 'include' and 'lib'.        
   
    Zheng Gong
    2020/05/06
