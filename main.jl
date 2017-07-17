# main.jl
# exact diagonalization study of J-K Hamiltonian
# Alan Morningstar
# May 2017


using ArgParse;
using DataFrames;


include("utils.jl");
include("lattice.jl");
include("basis.jl");
include("sparseHermitian.jl");
include("sparseHam.jl");
include("sparseS2.jl");


# main function
function main(Lx::Int64,Ly::Int64,J1::Float64,J2::Float64,K::Float64,Sz::Int64,mx::Int64,my::Int64)
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
    println("Sz = ",Sz);
    println("mx = ",mx);
    println("my = ",my);

    # number of eigenvalues desired
    numEigs::Int64 = 6;
    # a tolerance for error
    tolerance::Float64 = 10.^(-8.);
    # ritzVec = true if you want the eigenvectors returned too
    ritzVec::Bool = true;
    # number of Krylov vectors in eigenvalue calculation
    numKrylovVecs::Int64 = 40;
    # maximum number of iterations to converge eigenvalues
    maxIter::Int64 = 200;

    # plaquette (x,y) vectors, locating p1,p2,p3,p4 on the plaquette of the p1 site and p1D,p2D,p1L,p3L on adjacent plaquettes
    neighborVectors::Array{Tuple{Int64,Int64},1} = [(0,0),(1,0),(0,1),(1,1),(0,-1),(1,-1),(-1,0),(-1,1)];
    # define the lattice
    l::lattice = lattice(Lx,Ly,neighborVectors);

    # specify symmetry sector
    s::sector = sector(Sz,kx,ky);

    # construct the basis
    basis::SzkxkyBasis = SzkxkyBasis{UInt64}(l,s);
    println("Dimension of reduced Hilbert space is ",basis.dim,".");

    # couplings type to make passing J1,K easier
    c::couplings = couplings(J1,J2,K);

    # build the sparse Hamiltonian
    println("Building the Hamiltonian.");
    H = constructSparseHam(basis,c,s,l);

    # compute eigenvalues
    #:LM stands for largest magnitude, :SR for smallest real part
    println("Computing eigenvalues and eigenvectors.");
    eigsResult = eigs(H; nev=numEigs,ncv=numKrylovVecs,maxiter=maxIter, which=:SR, tol=tolerance, ritzvec=ritzVec);

    # clear Hamiltonian memory
    H = 0;

    # compile data
    println("Compiling data.");
    # energies
    EData::Array{Float64,1} = real(eigsResult[1]);
    # Sz values
    SzData::Array{Int64,1} = fill(Sz,numEigs);
    # mx values
    mxData::Array{Int64,1} = fill(mx,numEigs);
    # my values
    myData::Array{Int64,1} = fill(my,numEigs);
    # S(S+1) values
    S2Data::Array{Int64,1} = zeros(Int64,numEigs);
    # allocate memory before loop
    psi::Array{Complex128,1} = Array{Complex128,1}(basis.dim);
    S2psi::Array{Complex128,1} = Array{Complex128,1}(basis.dim);
    for i::Int64 in 1:numEigs
        psi = eigsResult[2][:,i];
        S2_mul_psi!(basis,s,l,S2psi,psi);
        S2Data[i] = round(Int64,real(dot(psi,S2psi))[1])
    end;

    # clear eigsResult memory
    eigsResult = 0;

    # create DataFrame
    df::DataFrame = DataFrame(E=EData,Ssqrd=S2Data,Sz=SzData,mx=mxData,my=myData);
    println(df);
    # sort by energy
    sort!(df);

    return df;

end;


#-- parse command line arguments
s = ArgParseSettings();

@add_arg_table s begin

    "--Lx"
        help = "length of lattice in x direction";
        arg_type = Int64;
        default = 4;
    "--Ly"
        help = "length of lattice in y direction";
        arg_type = Int64;
        default = 4;
    "--J1"
        help = "nearest neighbor coupling";
        arg_type = Float64;
        default = 1.0;
    "--J2"
        help = "next nearest neighbor coupling";
        arg_type = Float64;
        default = 0.0;
    "--K"
        help = "ring exchange coupling";
        arg_type = Float64;
        default = 0.0;
    "--Sz"
        help = "total z spin";
        arg_type = Int64;
        default = 0;
    "--mx"
        help = "x momentum mx such that kx = 2 pi mx / Lx";
        arg_type = Int64;
          default = 0;

    "--my"
        help = "y momentum my such that ky = 2 pi my / Ly";
        arg_type = Int64;
        default = 0;

end;

const argsDict = parse_args(s);


#-- execute main function
data = main(
            argsDict["Lx"],
            argsDict["Ly"],
            argsDict["J1"],
            argsDict["J2"],
            argsDict["K"],
            argsDict["Sz"],
            argsDict["mx"],
            argsDict["my"]
            );


#-- save data
dataFileName = "testing/Lx=" * string(argsDict["Lx"]) * "_Ly=" * string(argsDict["Ly"]) * "_J1=" * string(argsDict["J1"]) * "_J2=" * string(argsDict["J2"]) * "_K=" * string(argsDict["K"]) * "_Sz=" * string(argsDict["Sz"]) * "_mx=" * string(argsDict["mx"]) * "_my=" * string(argsDict["my"]) * ".csv";

writetable(dataFileName, data);
