# sparseS2.jl
# sparse representatin of S^2 operator
# Alan Morningstar
# June 2017


# function for building the sparse S^2 operator
function constructSparseS2(basis::SzkxkyBasis,s::sector,l::lattice)

    # lattice
    Lx::Int64 = l.Lx;
    Ly::Int64 = l.Ly;
    N::Int64 = l.N;
    # momentum
    kx::Float64 = s.kx;
    ky::Float64 = s.ky;

    # store location and value of non-zero matrix elements in CSC format
    Jpointers::Array{Int32,1} = Array{Int32,1}(basis.dim+1);
    I::Array{Int32,1} = Int32[];
    M::Array{Complex128,1} = Complex128[];

    # allocate memory once before the loops
    # -------------------------------------
    b::UInt64 = 0;
    sb::Array{UInt64,1} = Array{UInt64,1}(N);
    nb::Float64 = 0.0;
    S2bb::Complex128 = 0.0+0.0im;
    sPwij::Array{Int64,2} = Array{Int64,2}(N,N);
    a::UInt64 = UInt64(0);
    aRep::UInt64 = UInt64(0);
    lx::Int64 = ly::Int64 = 0;
    aRepIndex::Int32 = Int32(0);
    S2ij::Complex128 = 0.0+0.0im;
    # -------------------------------------

    # loop over basis
    for bIndex::Int64 in 1:basis.dim
        # CSC formatting
        Jpointers[bIndex] = Int32(length(I)+1);

        # spin states of basis
        b = basis.b[bIndex]; # integer rep.
        for bit in 1:l.N
            sb[bit] = readBit(b,bit); # spin array rep.
        end;
        # normalization constant
        nb = basis.n[bIndex];

        # fill the sPwij table
        for i::Int64 in 1:N
            for j::Int64 in 1:(i-1)
                sPwij[i,j] = simplePower(sb[i]+sb[j]);
            end;
        end;

        # initialize the diagonal matrix element
        S2bb = 0.75*N + 0.0im;

        # loop over pairs of lattice sites
        for i::Int64 in 1:N
            for j::Int64 in 1:(i-1)

                # contribute to the diagonal matrix element
                # -----------------------------------------
                S2bb += 0.5*sPwij[i,j];

                # compute off diagonal matrix elements
                # ------------------------------------

                if sb[i] != sb[j]
                    # the bra
                    a = XiXj(b,i,j);
                    # the rep and translation of the bra
                    aRep,lx,ly = representative(a,l);
                    # search for this rep in the basis
                    aRepIndex = basisIndex(aRep,basis);
                    if aRepIndex != 0
                        # the matrix element
                        S2ij = 0.5*(1-sPwij[i,j])*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);

                        push!(I,aRepIndex);
                        push!(M,S2ij);
                    end;
                end;

            end;
        end;

        # push diagonal matrix element to list of matrix elements
        push!(I,bIndex);
        push!(M,S2bb);

        # CSC formatting (turns out this is unnecessary)
        # sortTwo!(I[Jpointers[bIndex]:end],M[Jpointers[bIndex]:end],1,length(I[Jpointers[bIndex]:end]));

    end;

    # CSC formatting
    Jpointers[end] = Int32(length(I)+1);

    S2::SparseMatrixCSC{Complex128,Int32} = SparseMatrixCSC{Complex128,Int32}(basis.dim, basis.dim, Jpointers, I, M);
    return S2;
end;
