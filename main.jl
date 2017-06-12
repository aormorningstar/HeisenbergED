# main.jl
# exact diagonalization study of J-K Hamiltonian
# Alan Morningstar
# May 2017


include("utils.jl");
include("lattice.jl");
include("basis.jl");
include("sparseHam.jl")


# main function
function main(Lx::Int64,Ly::Int64,J1::Float64,K::Float64,n1::Int64,mx::Int64,my::Int64)
    # number of lattice sites
    N::Int64 = Lx*Ly;
    # momentum
    kx::Float64 = 2.0*pi*mx/Lx;
    ky::Float64 = 2.0*pi*my/Ly;

    # print details
    println("N = ",N);
    println("J = ",J1);
    println("K = ",K);
    println("n1 = ",n1);
    println("mx = ",mx);
    println("my = ",my);

    # number of eigenvalues desired
    numEigs::Int64 = 4;
    # a tolerance for error
    tolerance::Float64 = 10.^(-5.);
    # ritzVec = true if you want the eigenvectors returned too
    ritzVec::Bool = false;
    # number of Krylov vectors in eigenvalue calculation
    numKrylovVecs::Int64 = 20;
    # maximum number of iterations to converge eigenvalues
    maxIter::Int64 = 200;

    # plaquette (x,y) vectors, locating p1,p2,p3,p4 on the plaquette of the p1 site and p1D,p2D,p1L,p3L on adjacent plaquettes
    neighborVectors::Array{Tuple{Int64,Int64},1} = [(0,0),(1,0),(0,1),(1,1),(0,-1),(1,-1),(-1,0),(-1,1)];
    # define the lattice
    l::lattice = lattice(Lx,Ly,neighborVectors);

    # specify symmetry sector
    s::sector = sector(n1,kx,ky);

    # construct the basis
    basis::SzkxkyBasis = SzkxkyBasis(l,s);
    print("Dimension of reduced Hilbert space is ",basis.dim,". \n");

    # couplings type to make passing J1,K easier
    c::couplings = couplings(J1,K);

    # build the sparse Hamiltonian
    H::SparseMatrixCSC{Complex,Int64} = constructSparseHam(basis,c,s,l);

    # compute eigenvalues
    #:LM stands for largest magnitude, :SR for smallest real part
    eigsResult = eigs(H; nev=numEigs,ncv=numKrylovVecs,maxiter=maxIter, which=:SR, tol=tolerance, ritzvec=ritzVec);

    # print results
    print("Energies are: ",real(eigsResult[1]),". \n");

    # return energies
    return real(eigsResult[1]);
end;


#-- specify parameters of the calculation

# square lattice length
const Lx = 4;
const Ly = 4;
# NN coupling
const J1 = 1.0;
# plaquette coupling
const K = Float64[0.0];

# choose Sz sector by specifying number of 1s in basis states
const n1 = Int64[convert(Int64,(Lx*Ly)/2)];
# choose kx,ky by specifying mi such that mi is in 0:Li-1
const mx = Int64[0];
const my = Int64[0];


#-- execute main function

r1 = main(Lx,Ly,J1,K[1],n1[1],mx[1],my[1]);
