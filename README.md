# Exact Diagonalization of a Spin Hamiltonian
Arnoldi diagonalization of a spin Hamiltonian on a rectangular, two-dimensional lattice, implemented in julia. Uses Sz, translational, and spin-flip symmetry to reduce the Hilbert space dimension, and an ARPACK wrapper to find eigenvalues and eigenstates.

Outputs energies, and total-spin, spin-inversion, z-spin, and momentum quantum numbers for low lying eigenstates.

The provided spin Hamiltonian is the J1-J2-J3-J4-K model from Equation 4 of [Sachdev's paper (2017)](https://arxiv.org/abs/1705.06289), with J3 and J4 fixed to 0 currently.

For an almost-complete overview of methods used in this code, see Section 4.1 and 4.4 of [Sandvik's notes (2011)](https://arxiv.org/abs/1101.3281).

## Requirements
Julia 0.5, DataFrames, ArgParse

## Running a Simulation
There is a file named "main.jl" for running simulations. Command line arguments control Hamiltonian and lattice parameters, the symmetry sector of the Hamiltonian to be diagonalized, and the parameters of the Arnoldi diagonalization algorithm. Available command line arguments are:

| Argument        | Description           | Argument Type  |  Default Value |
| ------------- |-------------| -----| ----|
| Lx      | width of lattice | Int | 4 |
| Ly      | length of lattice      |   Int | 4 |
| J1 | n.n. coupling      |    Float | 1.0 |
| J2  | n.n.n. coupling  | Float  | 0.0  |
| K  | plaquette coupling  | Float  | 0.0  |
| Sz  | total spin along the z-axis  | Int  | 0  |
| mx  | momentum kx = 2 pi mx / Lx  | Int  | 0  |
| my  | momentum ky = 2 pi my / Ly  | Int  | 0  |
| z  | spin-inversion number  | Int  | 0  |
| verbose  | print details of simulation?  | Bool  | true  |
| numEigs  | number of eigenvalues calculated  | Int  | 4  |
| tol  | error tolerance of Arnoldi algorithm  | Float  | 10^(-8)  |
| numKrylovVecs  | dimension of the Krylov subspace  | Int  | 10  |
| maxIter  | maximum number of Arnoldi iterations  | Int  | 200  |

For example, a slurm job script could look like

```
#!/bin/bash                                                                        
#SBATCH --mem=100G                                                                                              
#SBATCH --time=2-00:00:00                                                                                       
julia main.jl --Lx 6 --Ly 6 --J1 1.0 --J2 0.0 --K 0.0 --Sz 0 --mx 0 --my 0 --z 1 --numEigs 4 --tol 0.00000001 --numKrylovVecs 10 --maxIter 200 --verbose true
```

where a 6x6 lattice is given 100 GB and 2 days to run (this should actually take ~85 GB and ~30 hours for the largest symmetry sectors of a 6x6 J-K Hamiltonian).

## Possible Improvements

*This section will be updated with more strategies for improvements to this code.

1. The current version uses julia's eigs(), a wrapper around ARPACK, which uses implicitely restarted Arnoldi iterations to find low lying eigenvalues. This method stores on the order of numKrylovVecs vectors, which is significantly larger than a Lanczos method, which only requires you to store three. So, need to implement a simple Lanczos algorithm such as that found in Section 4.2 of [(Sandvik, 2011)](https://arxiv.org/abs/1101.3281). This will also speed up the time to find eigenvalues because you can leave out any re-orthonormalization procedure and just make sure to monitor the Lanczos algorithm to avoid "ghost states" which arise when the Lanczos fails due to loss of orthogonality of the eigenvectors it is calculating.

## Author

* **Alan Morningstar**

