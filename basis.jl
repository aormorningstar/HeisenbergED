# basis.jl
# reduced basis in the Sz,kx,ky symmetry sector
# Alan Morningstar
# May 2017


# container for symmetry sector
immutable sector
    # number of 1 bits in the bit rep. of spin states in this sector
    n1::Int64;
    # momentum
    kx::Float64;
    ky::Float64;
end;


# function for computing the normalization constant Na such that the momentum state is 1/sqrt(Na) * sum(phase * TxTy |b>)
function normConstant(b::UInt64,l::lattice,s::sector)

    # sum of phases
    F::Complex128 = 0.0+0.0im;
    # normalization constant
    Na::Float64 = 0.0;
    # currrent translated state
    Tb::UInt64 = b;
    # set of translated states in integer rep.
    Tbs::Array{UInt64,1} = Array{UInt64}(l.N);

    # perform all translations
    for y in 0:l.Ly-1
        for x in 0:l.Lx-1
            # add to set of translated states
            Tbs[y*l.Lx+x+1] = Tb;

            if Tb == b
                # add to sum of phases
                F += exp(-1.0im*(s.kx*x+s.ky*y));
            elseif Tb < b
                return 0.0
            end;

            Tb = Tx(Tb,l);

        end;

        Tb = Ty(Tb,l);

    end;

    # compute and return normalization constant
    Na = abs2(F)*numUnique(Tbs);

    return Na;
end;


# the basis
immutable SzkxkyBasis
    # list of representatives of momentum basis states in integer representation
    b::Array{UInt64};
    # list of corresponding normalization constants
    n::Array{Float64};
    # dimension of Hilbert space
    dim::Int64;

    # constructor
    function SzkxkyBasis(l::lattice,s::sector)
        # initialize list of reps of momentum basis elements
        bList::Array{UInt64,1} = UInt64[];
        # initialize list of normalization constants
        nList::Array{Float64,1} = Float64[];

        # run up the binary odometer of Sz states
        # start with first n1 bits in state 1(down)
        b::UInt64 = 2^s.n1-1;

        while true
            # check if this state meets the required kx,ky
            n::Float64 = normConstant(b,l,s);

            # if valid rep state, add info to basis
            if n > 0.0
                push!(bList,b);
                push!(nList,n);
            end;

            # counter of 1 bits to the right (in binary convention)
            i::Int = 0;
            # position in bit array
            bit::Int = 1;
            while bit < l.N
                # find a 1 bit whose following neighbor is 0
                if readBit(b,bit) == 0x1
                    if readBit(b,bit+1) == 0x1
                        i += 1;
                    else
                        # shuffle bits over
                        for J in 1:i
                            b = setBit(b,J);
                        end;
                        for J in i+1:bit
                            b = clearBit(b,J);
                        end;
                        b = setBit(b,bit+1);
                        # then break to next loop
                        break
                    end;
                end;
                # to the next bit
                bit += 1;
            end;

            # if all 1s got shifted completely to the left, then all states are explored
            if bit == l.N
                break
            end;

        end;

        # construct basis
        new(bList,nList,length(bList));
    end;

end;


# search ordered basis for index of integer representation of spin state
function basisIndex(b::UInt64,basis::SzkxkyBasis)
    bIndex::UnitRange{Int64} = searchsorted(basis.b,b)::UnitRange{Int64};
    if !isempty(bIndex)
        # return Int32 because basis has less than 2 billion elements and need to save these in sparse Hamiltonian
        return Int32(bIndex[1])::Int32;
    else
        return Int32(0)::Int32;
    end;
end;
