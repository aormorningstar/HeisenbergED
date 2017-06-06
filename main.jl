# main.jl
# exact diagonalization study of J-K Hamiltonian
# Alan Morningstar
# May 2017


include("utils.jl");
include("lattice.jl");
include("basis.jl");
include("sparseHam.jl")


# main function
function main()
    #-- specify parameters of the calculation

    # square lattice length
    Lx::Int64 = 4;
    Ly::Int64 = 4;
    # number of sites
    N::Int64 = Lx*Ly;
    # NN coupling
    J1::Float64 = 1.0;
    # plaquette coupling
    K::Float64 = 0.0;

    # choose Sz sector by specifying number of 1s in basis states
    n1::Int64 = convert(Int,N/2);
    # choose kx,ky by specifying mi such that mi is in 0:Li-1
    mx::Int64 = 0;
    my::Int64 = 0;
    kx::Float64 = 2*pi*mx/Lx;
    ky::Float64 = 2*pi*my/Ly;

    # number of eigenvalues desired
    numEigs::Int64 = 4;
    # a tolerance for error
    tolerance::Float64 = 10.^(-5.);
    # ritzVec = true if you want the eigenvectors returned too
    ritzVec::Bool = false;
    # number of Krylov vectors in eigenvalue calculation
    numKrylovVecs::Int64 = 20;
    # maximum number of iterations to converge eigenvalues
    maxIter::Int64 = 300;

    # plaquette (x,y) vectors, locating p1,p2,p3,p4 on the plaquette of the p1 site and p1D,p2D,p1L,p3L on adjacent plaquettes
    neighborVectors::Array{Tuple{Int64,Int64},1} = [(0,0),(1,0),(0,1),(1,1),(0,-1),(1,-1),(-1,0),(-1,1)];
    # define the lattice
    l::lattice = lattice(Lx,Ly,neighborVectors);

    # specify symmetry sector
    s::sector = sector(n1,kx,ky);

    # construct the basis
    # 2.7 seconds for 6x4 lattice, basis is ~0.0002 GB
    basis::SzkxkyBasis = SzkxkyBasis(l,s);
    print("Dimension of reduced Hilbert space is ",basis.dim,". \n");

    # couplings type to make passing J1,K easier
    c::couplings = couplings(J1,K);

    # build the sparse Hamiltonian
    # takes 40 seconds for 6x4 J1-only Heisenberg model, H contains ~0.06 GB
    H::SparseMatrixCSC{Complex,Int64} = constructSparseHam(basis,c,s,l);

    # optional profiling
    # -----------------------------------
    # using ProfileView;

    # Profile.clear();
    # @profile constructSparseHam(basis,c,s,l);
    # ProfileView.view()
    # -----------------------------------

    #- compute eigenvalues
    # ~8.0 seconds for 6x4, numEigs=4, which=:SR, J1 Heisenberg Hamiltonian
    #:LM stands for largest magnitude, :SR for smallest real part
    eigsResult = eigs(H; nev=numEigs,ncv=numKrylovVecs,maxiter=maxIter, which=:SR, tol=tolerance, ritzvec=ritzVec);

    # print energies
    print("Energies are: ",real(eigsResult[1]),". \n");
end;

# execute main function
main();
