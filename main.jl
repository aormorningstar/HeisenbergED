# main.jl
# main function for ED study of J-K Hamiltonian
# Alan Morningstar
# June 2017


using ArgParse
using DataFrames


include("utils.jl")
include("lattice.jl")
include("basis.jl")
include("sparseHermitian.jl")
include("sparseHam.jl")
include("sparseS2.jl")


# main function
function main(
              Lx::Int64, # lattice length in x direction
              Ly::Int64, # lattice length in y direction
              J1::Float64, # nearest neighbor coupling
              J2::Float64, # next nearest neighbor coupling
              K::Float64, # plaquette coupling
              Sz::Int64, # spin along the z-axis
              mx::Int64, # momentum-related number, kx=2*pi*mx/Lx
              my::Int64, # momentum-related number, ky=2*pi*my/Ly
              z::Int64, # spin-inversion number
              numEigs::Int64, # number of eigenvalues desired
              tolerance::Float64, # strictness of convergence criterion
              numKrylovVecs::Int64, # number of vectors used in Krylov subspace
              maxIter::Int64, # max number of Arnoldi iterations
              verbose::Bool, # print details while running?
              )


    # number of lattice sites
    N = Lx*Ly
    # momentum
    kx = 2.0*pi*mx/Lx
    ky = 2.0*pi*my/Ly

    # print details
    verbose ?
        begin
            println("N = ",N)
            println("J1 = ",J1)
            println("J2 = ",J2)
            println("K = ",K)
            println("Sz = ",Sz)
            println("mx = ",mx)
            println("my = ",my)
            println("z = ",z)
        end : nothing

    # plaquette (x,y) vectors, locating p1,p2,p3,p4 on the plaquette of the p1 site and p1D,p2D,p1L,p3L on adjacent plaquettes
    neighborVectors = [(0,0),(1,0),(0,1),(1,1),(0,-1),(1,-1),(-1,0),(-1,1)]
    # define the lattice
    l = lattice(Lx,Ly,neighborVectors)

    # specify symmetry sector
    s = sector(Sz,kx,ky,z)

    # construct the basis
    verbose ? println("Constructing the basis.") : nothing
    basis = reducedBasis{UInt64}(l,s)
    verbose ? println("Dimension of reduced Hilbert space is ",basis.dim,".") : nothing

    # couplings type to make passing J1,J2,K easier
    c = couplings(J1,J2,K)

    # build the sparse Hamiltonian
    verbose ? println("Building the Hamiltonian.") : nothing
    H = constructSparseHam(basis,c,s,l)

    # compute eigenvalues
    # ritzVec = true if you want the eigenvectors returned too
    #:LM stands for largest magnitude, :SR for smallest real part
    verbose ? println("Computing eigenvalues and eigenvectors.") : nothing
    eigsResult = eigs(H; nev=numEigs,ncv=numKrylovVecs,maxiter=maxIter, which=:SR, tol=tolerance, ritzvec=true)

    # clear Hamiltonian memory manually
    H = nothing

    # compile data
    verbose ? println("Compiling data.") : nothing
    # energies
    EData = real(eigsResult[1])
    # Sz values
    SzData = fill(Sz,numEigs)
    # mx values
    mxData = fill(mx,numEigs)
    # my values
    myData = fill(my,numEigs)
    # z values
    zData = fill(z,numEigs);
    # S(S+1) values
    S2Data = round(Int64,real(S2expectations(basis,s,l,eigsResult[2])));

    # clear eigsResult memory manually
    eigsResult = nothing

    # create DataFrame
    df = DataFrame(E=EData,Ssqrd=S2Data,Sz=SzData,mx=mxData,my=myData,z=zData)
    verbose ? println(df) : nothing
    # sort by energy
    sort!(df)

    return df

end


#-- parse command line arguments
s = ArgParseSettings()

@add_arg_table s begin

    "--Lx"
        help = "length of lattice in x direction"
        arg_type = Int64
        default = 4
    "--Ly"
        help = "length of lattice in y direction"
        arg_type = Int64
        default = 4
    "--J1"
        help = "nearest neighbor coupling"
        arg_type = Float64
        default = 1.0
    "--J2"
        help = "next nearest neighbor coupling"
        arg_type = Float64
        default = 0.0
    "--K"
        help = "ring exchange coupling"
        arg_type = Float64
        default = 0.0
    "--Sz"
        help = "total z spin"
        arg_type = Int64
        default = 0
    "--mx"
        help = "x momentum mx such that kx = 2 pi mx / Lx"
        arg_type = Int64
        default = 0
    "--my"
        help = "y momentum my such that ky = 2 pi my / Ly"
        arg_type = Int64
        default = 0
    "--z"
        help = "z = +-1 spin-inversion quantum number"
        arg_type = Int64
        default = 1
    "--verbose"
        help = "print info during computation? true or false"
        arg_type = Bool
        default = true
    "--numEigs"
        help = "number of eigenvalues desired"
        arg_type = Int64
        default = 4
    "--tol"
        help = "tolerance for error in Arnoldi convergence"
        arg_type = Float64
        default = 10.^(-8.)
    "--numKrylovVecs"
        help = "number of vectors in the Krylov subspace"
        arg_type = Int64
        default = 10
    "--maxIter"
        help = "nmax number of Arnoldi iterations"
        arg_type = Int64
        default = 200

end

const argsDict = parse_args(s)


#-- execute main function
data = main(
            argsDict["Lx"],
            argsDict["Ly"],
            argsDict["J1"],
            argsDict["J2"],
            argsDict["K"],
            argsDict["Sz"],
            argsDict["mx"],
            argsDict["my"],
            argsDict["z"],
            argsDict["numEigs"],
            argsDict["tol"],
            argsDict["numKrylovVecs"],
            argsDict["maxIter"],
            argsDict["verbose"]
            )


#-- save data
dataFileName = "specData/Lx=" * string(argsDict["Lx"]) * "_Ly=" * string(argsDict["Ly"]) * "_J1=" * string(argsDict["J1"]) * "_J2=" * string(argsDict["J2"]) * "_K=" * string(argsDict["K"]) * "_Sz=" * string(argsDict["Sz"]) * "_mx=" * string(argsDict["mx"]) * "_my=" * string(argsDict["my"]) * "_z=" * string(argsDict["z"]) * ".csv"

writetable(dataFileName, data)
