# basis.jl
# reduced basis in the Sz,kx,ky,z symmetry sector
# Alan Morningstar
# May 2017


# container for symmetry sector
immutable sector

    # spin along the z-axis
    Sz::Int64
    # momentum in x and y directions
    kx::Float64
    ky::Float64
    # spin inversion quantum number
    z::Int64

end


# function for computing the normalization constant Na of
# the Tx,Ty,Z symmetrized state
function normConstant{I<:Integer}(b::I,l::lattice,s::sector)

    # sum of phases in the symmetrization
    F::Complex128 = 0.0+0.0im
    # currrent T and Z transformed states
    Tb::I = b
    ZTb::I = Z(b,l)
    # set of transformed states in integer rep.
    Tbs::Vector{I} = Vector{I}()

    # perform all transformations
    for y::Int64 in 0:l.Ly-1
        for x::Int64 in 0:l.Lx-1

            # add to set of translated states
            appendSet!(Tb,Tbs)

            if Tb == b
                # add to sum of phases
                F += exp(-1.0im*(s.kx*x+s.ky*y))
            elseif Tb < b
                return 0.0
            end

            # spin-invert the state
            ZTb = Z(Tb,l)

            # add inverted state to set of states
            appendSet!(ZTb,Tbs)

            if ZTb == b
                # add to sum of phases
                F += s.z*exp(-1.0im*(s.kx*x+s.ky*y))
            elseif ZTb < b
                return 0.0
            end

            # translate state in x direction
            Tb = Tx(Tb,l)

        end

        # translate state in y direction
        Tb = Ty(Tb,l)

    end

    # compute and return normalization constant
    return abs2(F)*length(Tbs)

end


# the basis
immutable reducedBasis{I<:Integer}

    # list of representatives of momentum basis states in integer representation
    b::Array{I,1}
    # list of corresponding normalization constants
    n::Array{Float64,1}
    # dimension of reduced Hilbert space
    dim::Int64

    # constructor
    reducedBasis(l::lattice,s::sector) =
    begin
        # initialize list of reps of momentum basis elements
        bList::Array{I,1} = Array{I,1}()
        # initialize list of normalization constants
        nList::Array{Float64,1} = Float64[]

        # run up the binary odometer of Sz states
        # start with first n1 bits in state 1(down)
        n1::Int64 = div(l.N,2)-s.Sz
        b::I = convert(I,2^n1-1)

        # allocate memory once here
        n::Float64 = 0.0
        i::Int64 = 0
        bit::Int64 = 0

        while true
            # check if this state has normalization constant = 0, if yes
            # then add it to the basis
            n = normConstant(b,l,s)

            if n > 0.000001
                push!(bList,b)
                push!(nList,n)
            end

            # counter of 1 bits to the right (in binary convention)
            # NOTE: these are just temporary variables for the 'binary odometer'
            i = 0
            # position in bit array
            bit = 1

            while bit < l.N

                # find a 1 bit whose following neighbor is 0
                if readBit(b,bit) == 1

                    if readBit(b,bit+1) == 1
                        i += 1
                    else

                        # shuffle bits over
                        for J::Int64 in 1:i
                            b = setBit(b,J)
                        end

                        for J::Int64 in i+1:bit
                            b = clearBit(b,J)
                        end

                        b = setBit(b,bit+1)

                        # then break to next loop
                        break

                    end

                end

                # to the next bit
                bit += 1

            end

            # if all 1s got shifted completely to the left, then all states
            # are explored
            if bit == l.N
                break
            end

        end

        # construct basis
        new(bList,nList,length(bList))

    end

end


# search ordered basis for index (or location in the basis) of an integer
# representation of a spin state
# NOTE: basis has less than ~ 2 billion states, else modify to use Int64
function basisIndex{I<:Integer}(b::I,basis::reducedBasis{I})

    bIndex::UnitRange{Int64} = searchsorted(basis.b::Array{I,1},b)::UnitRange{Int64}

    if !isempty(bIndex)
        # return Int32 because basis has less than 2 billion elements
        # and need to save these in sparse Hamiltonian
        return Int32(bIndex[1])::Int32
    else
        return Int32(0)::Int32
    end

end


# return the type of the basis elements
Base.eltype(basis::reducedBasis) = eltype(basis.b)
