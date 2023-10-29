# gpu-diagonalizer
This repository contains a modified version of the original diag_lapack.f library designed to exploit GPUs for performing linear algebra operations, including matrix diagonalization, matrix-matrix, and matrix-vector multiplications within the CRYSTAL code. The offloading strategy reported here is based on OpenACC directives. This file is currently under development and has only been tested on NVIDIA Ampere GPUs. It is intended to replace the file with the same name in the public version of the CRYSTAL23 code.

# Leonardo GPU partition scaling
![Alt text](Leonardo_timing_SiO2_X4.png "Leonardo GPU partition scaling")
*Leonardo booster module computing time using one and two nodes, for a supercell of SiO2 containing 648 atoms and 13392 atomic orbitals. SCF performed at PBE/m-6-311G(d)_Heyd_2005 level of theory.*
