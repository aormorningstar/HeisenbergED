# Exact Diagonalization of a Spin Hamiltonian
Arnoldi diagonalization of a spin Hamiltonian on a rectangular, two-dimensional lattice, implemented in julia. Uses Sz, translational, and spin-flip symmetry to reduce the Hilbert space dimension, and an ARPACK wrapper to find eigenvalues and eigenstates.

Outputs energies, and total-spin, spin-inversion, z-spin, and momentum quantum numbers for low lying eigenstates.

The provided spin Hamiltonian is the J1-J2-J3-J4-K model from Equation 4 of https://arxiv.org/abs/1705.06289, with J3 and J4 fixed to 0 currently.

## Requirements
Julia 0.5, DataFrames, ArgParse

## Running a Simulation
There is a file named "main.jl" for running simulations. Command line arguments control Hamiltonian and lattice parameters, the symmetry sector of the Hamiltonian to be diagonalized, and the parameters of the Arnoldi diagonalization algorithm. Available command line arguments are:

--Lx\\
    help = "length of lattice in x direction"\\
    arg_type = Int64\\
    default = 4\\

--Ly\\
    help = "length of lattice in y direction"\\
    arg_type = Int64\\
    default = 4\\

--J1\\
    help = "nearest neighbor coupling"\\
    arg_type = Float64\\
    default = 1.0\\

--J2\\
    help = "next nearest neighbor coupling"\\
    arg_type = Float64\\
    default = 0.0\\

--K\\
    help = "ring exchange coupling"\\
    arg_type = Float64\\
    default = 0.0\\

--Sz\\
    help = "total z spin"\\
    arg_type = Int64\\
    default = 0\\

--mx\\
    help = "x momentum mx such that kx = 2 pi mx / Lx"\\
    arg_type = Int64\\
    default = 0\\

--my\\
    help = "y momentum my such that ky = 2 pi my / Ly"\\
    arg_type = Int64\\
    default = 0\\

--z\\
    help = "z = +-1 spin-inversion quantum number"\\
    arg_type = Int64\\
    default = 1\\

--verbose\\
    help = "print info during computation? true or false"\\
    arg_type = Bool\\
    default = true\\

--numEigs
    help = "number of eigenvalues desired"\\
    arg_type = Int64\\
    default = 4\\

--tol\\
    help = "tolerance for error in Arnoldi convergence"\\
    arg_type = Float64\\
    default = 10.^(-8.)\\

--numKrylovVecs\\
    help = "number of vectors in the Krylov subspace"\\
    arg_type = Int64\\
    default = 10\\

--maxIter\\
    help = "nmax number of Arnoldi iterations"\\
    arg_type = Int64\\
    default = 200.\\

For example, a slurm job script could like like

\#!/bin/bash\\                                                                        
\#SBATCH --mem=100G\\                                                                                                
\#SBATCH --time=2-00:00:00\\                                                                                            
julia main.jl --Lx 6 --Ly 6 --J1 1.0 --J2 0.0 --K 0.0 --Sz 0 --mx 0 --my 0 --z 1 --numEigs 4 --tol 0.00000001 --numKrylovVecs 10 --maxIter 200 --verbose true

where a 6x6 lattice is givven 100 GB and 2 days to run (this should actually take ~85 GB and ~30 hours for the largest symmetry sectors of a 6x6 J-K Hamiltonian).


### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc

