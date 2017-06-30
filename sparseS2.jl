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

    # basis type
    bType::Type = eltype(basis);

    # store location and value of non-zero matrix elements in CSC format
    Jpointers::Vector{Int32} = Array{Int32}(basis.dim+1);
    I::Vector{Int32} = Int32[];
    M::Vector{Complex128} = Complex128[];
    Mpointers::Vector{Int32} = Int32[];

    # allocate memory once before the loops
    # -------------------------------------
    b::bType = 0;
    sb::Array{bType,1} = Array{bType,1}(N);
    nb::Float64 = 0.0;
    S2bb::Complex128 = 0.0+0.0im;
    sPwij::Array{Int64,2} = Array{Int64,2}(N,N);
    a::bType = convert(bType,0);
    aRep::bType = convert(bType,0);
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
        # NOTE: the leading factor of 0.5 in each diagonal term is due to the Hermitian sparse matrix which requires diagonal elements to be divided by 2 before being passed
        S2bb = 0.5*(0.75*N + 0.0im); # see NOTE

        # loop over pairs of lattice sites
        for i::Int64 in 1:N
            for j::Int64 in 1:(i-1)

                # contribute to the diagonal matrix element
                # -----------------------------------------
                S2bb += 0.5*0.5*sPwij[i,j]; # see NOTE

                # compute off diagonal matrix elements
                # ------------------------------------

                if sb[i] != sb[j]
                    # the bra
                    a = XiXj(b,i,j);
                    # the rep and translation of the bra
                    aRep,lx,ly = representative(a,l);
                    # search for this rep in the basis
                    aRepIndex = basisIndex(aRep,basis);
                    if aRepIndex != 0 && bIndex > aRepIndex # only keep upper triangle
                        # the matrix element
                        S2ij = 0.5*(1-sPwij[i,j])*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);

                        push!(I,aRepIndex);
                        # push!(M,S2ij);
                        push!(Mpointers,appendSet!(S2ij,M));
                    end;
                end;

            end;
        end;

        # push diagonal matrix element to list of matrix elements
        push!(I,bIndex);
        # push!(M,S2bb);
        push!(Mpointers,appendSet!(S2bb,M));

        # CSC formatting (turns out this is unnecessary)
        # sortTwo!(I[Jpointers[bIndex]:end],M[Jpointers[bIndex]:end],1,length(I[Jpointers[bIndex]:end]));

    end;

    # CSC formatting
    Jpointers[end] = Int32(length(I)+1);

    # S2::SparseMatrixCSC{Complex128,Int32} = SparseMatrixCSC{Complex128,Int32}(basis.dim, basis.dim, Jpointers, I, M);
    # S2::sparseHermitian{Int32,Complex128} = sparseHermitian{Int32,Complex128}(basis.dim,basis.dim,Jpointers,I,M);
    S2::sparseHermitian{Int32,Complex128} = sparseHermitian{Int32,Complex128}(basis.dim,basis.dim,Jpointers,I,M,Mpointers);

    return S2;
end;


# function for in place multiplication by the S^2 operator
function S2_mul_psi!(basis::SzkxkyBasis,s::sector,l::lattice,S2psi::Vector{Complex128},psi::Vector{Complex128})

    # clear output Vector
    S2psi .= 0.0+0.0im;

    # lattice
    Lx::Int64 = l.Lx;
    Ly::Int64 = l.Ly;
    N::Int64 = l.N;
    # momentum
    kx::Float64 = s.kx;
    ky::Float64 = s.ky;

    # basis type
    bType::Type = eltype(basis);

    # allocate memory once before the loops
    # -------------------------------------
    b::bType = 0;
    sb::Array{bType,1} = Array{bType,1}(N);
    nb::Float64 = 0.0;
    S2bb::Complex128 = 0.0+0.0im;
    sPwij::Array{Int64,2} = Array{Int64,2}(N,N);
    a::bType = convert(bType,0);
    aRep::bType = convert(bType,0);
    lx::Int64 = ly::Int64 = 0;
    aRepIndex::Int32 = Int32(0);
    S2ij::Complex128 = 0.0+0.0im;
    # -------------------------------------

    # loop over basis
    for bIndex::Int64 in 1:basis.dim

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
                        # contributing to the output vector
                        S2psi[aRepIndex] += S2ij*psi[bIndex];
                    end;
                end;

            end;
        end;

        # diagonal element contribution to output
        S2psi[bIndex] += S2bb*psi[bIndex];

    end;

    return S2psi::Vector{Complex128};
end;
