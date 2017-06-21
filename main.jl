# main.jl
# exact diagonalization study of J-K Hamiltonian
# Alan Morningstar
# May 2017

using DataFrames;

include("utils.jl");
include("lattice.jl");
include("basis.jl");
include("sparseHam.jl");
# include("sparseS2.jl");


# main function
function main(Lx::Int64,Ly::Int64,J1::Float64,J2::Float64,K::Float64,n1::Int64,mx::Int64,my::Int64)
    # number of lattice sites
    N::Int64 = Lx*Ly;
    # momentum
    kx::Float64 = 2.0*pi*mx/Lx;
    ky::Float64 = 2.0*pi*my/Ly;

    # print details
    println("N = ",N);
    println("J1 = ",J1);
    println("J2 = ",J2);
    println("K = ",K);
    println("n1 = ",n1);
    println("mx = ",mx);
    println("my = ",my);

    # number of eigenvalues desired
    numEigs::Int64 = 6;
    # a tolerance for error
    tolerance::Float64 = 10.^(-5.);
    # ritzVec = true if you want the eigenvectors returned too
    ritzVec::Bool = true;
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
    println("Dimension of reduced Hilbert space is ",basis.dim,".");

    # couplings type to make passing J1,K easier
    c::couplings = couplings(J1,J2,K);

    # build the sparse Hamiltonian
    println("Building the Hamiltonian.");
    H::SparseMatrixCSC{Complex,Int32} = constructSparseHam(basis,c,s,l);

    # compute eigenvalues
    #:LM stands for largest magnitude, :SR for smallest real part
    println("Computing eigenvalues and eigenvectors.");
    eigsResult = eigs(H; nev=numEigs,ncv=numKrylovVecs,maxiter=maxIter, which=:SR, tol=tolerance, ritzvec=ritzVec);

    # build the sparse S^2 operator
    # println("Building the S^2 operator.");
    # S2::SparseMatrixCSC{Complex,Int32} = constructSparseS2(basis,s,l);

    # compile data
    println("Compiling data.");
    # energies
    EData::Array{Float64,1} = real(eigsResult[1]);
    # Sz values
    SzData::Array{Float64,1} = fill((N-2*n1)/2,numEigs);
    # mx values
    mxData::Array{Int64,1} = fill(mx,numEigs);
    # my values
    myData::Array{Int64,1} = fill(my,numEigs);
    # S(S+1) values
    # S2Data::Array{Int64,1} = zeros(Int64,numEigs);
    # println("Compiling data 2.");
    # for i::Int64 in 1:numEigs
    #   psi::Array{Complex128,1} = eigsResult[2][:,i];
    #   println("Compiling data 3.");
    #   S2Data[i] = round(Int64,real(dot(psi,S2*psi))[1])
    # end;

    # create DataFrame
    println("Filling DataFrame.");
    df::DataFrame = DataFrame(E=EData,Sz=SzData,mx=mxData,my=myData);
    #,Ssqrd=S2Data);
    println(df);

    return df;

end;


#-- specify parameters of the calculation

# square lattice length
const Lx = 6;
const Ly = 4;
# NN coupling
const J1 = 1.0;
# NNN coupling
const J2 = 0.0;
# plaquette coupling
const K = 0.0;

# choose Sz sector by specifying number of 1s in basis states
const n1List = convert(Int64,(Lx*Ly)/2):(convert(Int64,(Lx*Ly)/2)+8);
# choose kx,ky by specifying mi such that mi is in 0:Li-1
const mxList = 0:(Lx-1);
const myList = 0:(Ly-1);


#-- execute main function for all parameter values and collect data
data = DataFrame();

for n1 in n1List
    for mx in mxList
        for my in myList
            df = main(Lx,Ly,J1,J2,K,n1[1],mx[1],my[1]);
            data = vcat(data,df);
        end;
    end;
end;

# write data
dataFileName = "specData/Lx=" * string(Lx) * "_Ly=" * string(Ly) * "_J1=" * string(J1) * "_J2=" * string(J2) * "_K=" * string(K) * ".csv";
writetable(dataFileName, data);
