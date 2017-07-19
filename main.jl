# main.jl
# exact diagonalization study of J-K Hamiltonian
# Alan Morningstar
# May 2017


using ArgParse
using DataFrames


include("utils.jl")
include("lattice.jl")
include("basis.jl")
include("sparseHermitian.jl")
include("sparseHam.jl")
include("sparseS2.jl")


# main function
function main(Lx::Int64,Ly::Int64,J1::Float64,J2::Float64,K::Float64,Sz::Int64,mx::Int64,my::Int64,z::Int64,verbose::Bool)
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

    # number of eigenvalues desired
    numEigs = 4
    # a tolerance for error
    tolerance = 10.^(-8.)
    # ritzVec = true if you want the eigenvectors returned too
    ritzVec = true
    # number of Krylov vectors in eigenvalue calculation
    numKrylovVecs = 40
    # maximum number of iterations to converge eigenvalues
    maxIter = 200

    # plaquette (x,y) vectors, locating p1,p2,p3,p4 on the plaquette of the p1 site and p1D,p2D,p1L,p3L on adjacent plaquettes
    neighborVectors = [(0,0),(1,0),(0,1),(1,1),(0,-1),(1,-1),(-1,0),(-1,1)]
    # define the lattice
    l = lattice(Lx,Ly,neighborVectors)

    # specify symmetry sector
    s = sector(Sz,kx,ky,z)

    # construct the basis
    basis = reducedBasis{UInt64}(l,s)
    verbose ? println("Dimension of reduced Hilbert space is ",basis.dim,".") : nothing

    # couplings type to make passing J1,J2,K easier
    c = couplings(J1,J2,K)

    # build the sparse Hamiltonian
    verbose ? println("Building the Hamiltonian.") : nothing
    H = constructSparseHam(basis,c,s,l)

    # compute eigenvalues
    #:LM stands for largest magnitude, :SR for smallest real part
    verbose ? println("Computing eigenvalues and eigenvectors.") : nothing
    eigsResult = eigs(H; nev=numEigs,ncv=numKrylovVecs,maxiter=maxIter, which=:SR, tol=tolerance, ritzvec=ritzVec)

    # clear Hamiltonian memory
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

    # clear eigsResult memory
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
            argsDict["verbose"]
            )


#-- save data
dataFileName = "testing/Lx=" * string(argsDict["Lx"]) * "_Ly=" * string(argsDict["Ly"]) * "_J1=" * string(argsDict["J1"]) * "_J2=" * string(argsDict["J2"]) * "_K=" * string(argsDict["K"]) * "_Sz=" * string(argsDict["Sz"]) * "_mx=" * string(argsDict["mx"]) * "_my=" * string(argsDict["my"]) * "_z=" * string(argsDict["z"]) * ".csv"

writetable(dataFileName, data)
