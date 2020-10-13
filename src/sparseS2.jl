# sparseS2.jl
# sparse S^2 operator
# Alan Morningstar
# June 2017


# function which takes in a matrix E whose columns are vectors and returns a
# list of S^2 expectation values for those column vectors
function S2expectations(basis::reducedBasis,s::sector,l::lattice,E::Array{Complex128,2})

    # lattice
    Lx::Int64 = l.Lx
    Ly::Int64 = l.Ly
    N::Int64 = l.N
    # momentum
    kx::Float64 = s.kx
    ky::Float64 = s.ky
    # spin-inversion number
    z::Float64 = s.z

    # list of S^2 expectations to be filled
    numVecs::Int64 = size(E)[2]
    S2s = zeros(Complex128,numVecs)

    # basis type
    bType::Type = eltype(basis)

    # allocate memory once before the loops
    # -------------------------------------
    b::bType = 0
    sb::Array{bType,1} = Array{bType,1}(N)
    nb::Float64 = 0.0
    S2bb::Complex128 = 0.0+0.0im
    sPwij::Array{Int64,2} = Array{Int64,2}(N,N)
    a::bType = convert(bType,0)
    aRep::bType = convert(bType,0)
    lx::Int64 = ly::Int64 = g::Int64 = 0
    aRepIndex::Int32 = Int32(0)
    S2ij::Complex128 = 0.0+0.0im
    # -------------------------------------

    # loop over basis
    for bIndex::Int64 in 1:basis.dim

        # spin states of basis
        b = basis.b[bIndex] # integer rep.
        for bit in 1:l.N
            sb[bit] = readBit(b,bit) # spin array rep.
        end

        # normalization constant
        nb = basis.n[bIndex]

        # fill the sPwij table
        for i::Int64 in 1:N
            for j::Int64 in 1:(i-1)
                sPwij[i,j] = simplePower(sb[i]+sb[j])
            end
        end

        # initialize the diagonal matrix element
        S2bb = 0.75*N + 0.0im

        # loop over pairs of lattice sites
        for i::Int64 in 1:N
            for j::Int64 in 1:(i-1)

                # contribute to the diagonal matrix element
                # -----------------------------------------
                S2bb += 0.5*sPwij[i,j]

                # compute off diagonal matrix elements
                # ------------------------------------

                if sb[i] != sb[j]
                    # the bra
                    a = XiXj(b,i,j)
                    # the rep and translation of the bra
                    aRep,lx,ly,g = representative(a,l)
                    # search for this rep in the basis
                    aRepIndex = basisIndex(aRep,basis)
                    if aRepIndex != 0
                        # the matrix element
                        S2ij = 0.5*(1-sPwij[i,j])*(z^g)*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb)
                        # contributing to the output
                        for vec in 1:numVecs
                            S2s[vec] += E[bIndex,vec]*conj(E[aRepIndex,vec])*S2ij
                        end
                    end
                end

            end
        end

        # contributing to the output from the diagonal matrix elements
        for vec in 1:numVecs
            S2s[vec] += E[bIndex,vec]*conj(E[bIndex,vec])*S2bb
        end

    end

    return S2s::Vector{Complex128}
end
