# sparseHam.jl
# sparse representatin of J-K Hamiltonian
# Alan Morningstar
# May 2017


# container for Hamiltonian couplings
immutable couplings
    # 1st nearest neighbor
    J1::Float64;
    # 2nd nearest neighbor
    J2::Float64;
    # plaquette terms
    K::Float64;
end;

# function for building the sparse Hamiltonian
function constructSparseHam(basis::SzkxkyBasis,c::couplings,s::sector,l::lattice)

    # lattice
    Lx::Int64 = l.Lx;
    Ly::Int64 = l.Ly;
    N::Int64 = l.N;
    # couplings
    J1::Float64 = c.J1;
    J2::Float64 = c.J2;
    K::Float64 = c.K;
    # momentum
    kx::Float64 = s.kx;
    ky::Float64 = s.ky;

    # non-zero couplings
    J1nz::Bool = (J1 != 0.0);
    J2nz::Bool = (J2 != 0.0);
    Knz::Bool = (K != 0.0);

    # basis type
    bType::Type = eltype(basis);

    # store location and value of non-zero matrix elements in CSC format
    Jpointers::Vector{Int32} = Vector{Int32}(basis.dim+1);
    I::Vector{Int32} = Int32[];
    M::Vector{Complex128} = Complex128[];
    Mpointers::Vector{Int32} = Int32[];

    # allocate memory once before the loops
    # -------------------------------------
    # nnz::Int64 = 0;
    b::bType = 0;
    sb::Array{bType,1} = Array{bType,1}(l.N);
    nb::Float64 = 0.0;
    Hbb::Complex128 = 0.0+0.0im;
    p::Array{Int64,1} = Array{Int64,1}(length(l.nbrs[1]));
    sp::Array{bType,1} = Array{UInt64,1}(length(p));
    sPw12::Int64 = sPw34::Int64 = sPw13::Int64 = sPw24::Int64 = sPw14::Int64 = sPw23::Int64 = sPw1234::Int64 = sPw1D2D::Int64 = sPw1L3L::Int64 = 0;
    c12::Bool = c13::Bool = c14::Bool = c23::Bool = c1234::Bool = true;
    a::bType = convert(bType,0);
    aRep::bType = convert(bType,0);
    lx::Int64 = ly::Int64 = 0;
    aRepIndex::Int32 = Int32(0);
    Hsite::Complex128 = 0.0+0.0im;
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

        # initialize the diagonal matrix element
        Hbb = 0.0;

        # loop over lattice sites
        for site::Int64 in 1:N

            # spins at nearby sites p1,p2,p3,p4,p1D,p2D,p1L,p3L
            p = l.nbrs[site];
            sp = sb[p];
            # some common factors in the matrix elements
            sPw12 = simplePower(sp[1] + sp[2]);
            sPw34 = simplePower(sp[3] + sp[4]);
            sPw13 = simplePower(sp[1] + sp[3]);
            sPw24 = simplePower(sp[2] + sp[4]);
            sPw14 = simplePower(sp[1] + sp[4]);
            sPw23 = simplePower(sp[2] + sp[3]);
            sPw1234 = simplePower(sp[1] + sp[2] + sp[3] + sp[4]);
            sPw1D2D = simplePower(sp[5] + sp[6]);
            sPw1L3L = simplePower(sp[7] + sp[8]);

            # some conditions
            c12 = ((sPw12 == -1) && (Knz || J1nz));
            c13 = ((sPw13 == -1) && (Knz || J1nz));
            c14 = ((sPw14 == -1) && (Knz || J2nz));
            c23 = ((sPw23 == -1) && (Knz || J2nz));
            c1234 = ((sp[1] + sp[2] + sp[3] + sp[4] == 2) && Knz);


            # contribute to the diagonal matrix element
            # -----------------------------------------
            if J1nz
                Hbb += 0.25*J1*( sPw12 + sPw13 );
            end;
            if J2nz
                Hbb += 0.25*J2*( sPw14 + sPw23 );
            end;
            if Knz
                Hbb += 0.125*K*sPw1234;
            end;

            # compute off diagonal matrix elements
            # ------------------------------------

            # the 12 term
            if c12
                # the bra
                a = XiXj(b,p[1],p[2]);
                # the rep and transaltion of the bra
                aRep,lx,ly = representative(a,l);
                # search for this rep in the basis
                aRepIndex = basisIndex(aRep,basis);
                if aRepIndex != 0 && (bIndex > aRepIndex) # only keep upper triangle
                    # the matrix element
                    Hsite = (0.5*J1+0.125*K*(sPw34-sPw1234+2*sPw1D2D))*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);

                    push!(I,aRepIndex);
                    push!(Mpointers,appendSet!(Hsite,M));
                elseif bIndex == aRepIndex
                    Hbb += 0.5*J1+0.125*K*(sPw34-sPw1234+2*sPw1D2D);
                end;
            end;

            # the 13 term
            if c13
                # the bra
                a = XiXj(b,p[1],p[3]);
                # the rep and transaltion of the bra
                aRep,lx,ly = representative(a,l);
                # search for this rep in the basis
                aRepIndex = basisIndex(aRep,basis);
                if aRepIndex != 0 && (bIndex > aRepIndex) # only keep upper triangle
                    # the matrix element
                    Hsite = (0.5*J1+0.125*K*(sPw24-sPw1234+2*sPw1L3L))*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);

                    push!(I,aRepIndex);
                    push!(Mpointers,appendSet!(Hsite,M));
                elseif bIndex == aRepIndex
                    Hbb += 0.5*J1+0.125*K*(sPw24-sPw1234+2*sPw1L3L);
                end;
            end;

            # the 14 term
            if c14
                # the bra
                a = XiXj(b,p[1],p[4]);
                # the rep and transaltion of the bra
                aRep,lx,ly = representative(a,l);
                # search for this rep in the basis
                aRepIndex = basisIndex(aRep,basis);
                if aRepIndex !=0 && (bIndex > aRepIndex) # only keep upper triangle
                    # the matrix element
                    Hsite = (0.5*J2-0.25*K*sPw23)*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);

                    push!(I,aRepIndex);
                    push!(Mpointers,appendSet!(Hsite,M));
                elseif bIndex == aRepIndex
                    Hbb += 0.5*J2-0.25*K*sPw23;
                end;
            end;

            # the 23 term
            if c23
                # the bra
                a = XiXj(b,p[2],p[3]);
                # the rep and transaltion of the bra
                aRep,lx,ly = representative(a,l);
                # search for this rep in the basis
                aRepIndex = basisIndex(aRep,basis);
                if aRepIndex != 0 && (bIndex > aRepIndex) # only keep upper triangle
                    # the matrix element
                    Hsite= (0.5*J2-0.25*K*sPw14)*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);

                    push!(I,aRepIndex);
                    push!(Mpointers,appendSet!(Hsite,M));
                elseif bIndex == aRepIndex
                    Hbb += 0.5*J2-0.25*K*sPw14;
                end;
            end;

            # the 1234 term
            if c1234
                # the bra
                a = XiXjXkXl(b,p[1],p[2],p[3],p[4]);
                # the rep and transaltion of the bra
                aRep,lx,ly = representative(a,l);
                # search for this rep in the basis
                aRepIndex = basisIndex(aRep,basis);
                if aRepIndex != 0 && (bIndex > aRepIndex) # only keep upper triangle
                    # the matrix element
                    Hsite = 0.125*K*(2-sPw12-sPw34-sPw13-sPw24+sPw14+sPw23)*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);

                    push!(I,aRepIndex);
                    push!(Mpointers,appendSet!(Hsite,M));
                elseif bIndex == aRepIndex
                    Hbb += 0.125*K*(2-sPw12-sPw34-sPw13-sPw24+sPw14+sPw23)
                end;
            end;

        end;

        # push diagonal matrix element to list of matrix elements
        # NOTE: the leading factor of 0.5 in each diagonal term is due to the Hermitian sparse matrix which requires diagonal elements to be divided by 2 before being passed
        push!(I,bIndex);
        push!(Mpointers,appendSet!(0.5*Hbb,M)); # see NOTE

    end;

    # CSC formatting
    Jpointers[end] = Int32(length(I)+1);

    H::sparseHermitian{Int32,Complex128} = sparseHermitian{Int32,Complex128}(basis.dim,basis.dim,Jpointers,I,M,Mpointers);
    return H;
end;


## THE REST OF THE FILE IS FOR TESTING ##

# # function for building the sparse Hamiltonian
# function constructSparseHam(basis::SzkxkyBasis,c::couplings,s::sector,l::lattice)
#
#     # lattice
#     Lx::Int64 = l.Lx;
#     Ly::Int64 = l.Ly;
#     N::Int64 = l.N;
#     # couplings
#     J1::Float64 = c.J1;
#     J2::Float64 = c.J2;
#     K::Float64 = c.K;
#     # momentum
#     kx::Float64 = s.kx;
#     ky::Float64 = s.ky;
#
#     # non-zero couplings
#     J1nz::Bool = (J1 != 0.0);
#     J2nz::Bool = (J2 != 0.0);
#     Knz::Bool = (K != 0.0);
#
#     # basis type
#     bType::Type = eltype(basis);
#
#     # store location and value of non-zero matrix elements in CSC format
#     Jpointers::Vector{Int32} = Vector{Int32}(basis.dim+1);
#     I::Vector{Int32} = Int32[];
#     M::Vector{Complex128} = Complex128[];
#     Mpointers::Vector{Int32} = Int32[];
#
#     # allocate memory once before the loops
#     # -------------------------------------
#     nnz::Int64 = 0;
#     b::bType = 0;
#     sb::Array{bType,1} = Array{bType,1}(l.N);
#     nb::Float64 = 0.0;
#     Hbb::Complex128 = 0.0+0.0im;
#     p::Array{Int64,1} = Array{Int64,1}(length(l.nbrs[1]));
#     sp::Array{bType,1} = Array{UInt64,1}(length(p));
#     sPw12::Int64 = sPw34::Int64 = sPw13::Int64 = sPw24::Int64 = sPw14::Int64 = sPw23::Int64 = sPw1234::Int64 = sPw1D2D::Int64 = sPw1L3L::Int64 = 0;
#     c12::Bool = c13::Bool = c14::Bool = c23::Bool = c1234::Bool = true;
#     a::bType = convert(bType,0);
#     aRep::bType = convert(bType,0);
#     lx::Int64 = ly::Int64 = 0;
#     aRepIndex::Int32 = Int32(0);
#     Hsite::Complex128 = 0.0+0.0im;
#     # -------------------------------------
#
#     # loop over basis
#     for bIndex::Int64 in 1:basis.dim
#         # CSC formatting
#         Jpointers[bIndex] = Int32(length(I)+1);
#
#         # count non-zero values
#         nnz = 0;
#
#         # spin states of basis
#         b = basis.b[bIndex]; # integer rep.
#         for bit in 1:l.N
#             sb[bit] = readBit(b,bit); # spin array rep.
#         end;
#         # normalization constant
#         nb = basis.n[bIndex];
#
#         # initialize the diagonal matrix element
#         Hbb = 0.0 + 0.0im;
#
#         # loop over lattice sites
#         for site::Int64 in 1:N
#
#             # spins at nearby sites p1,p2,p3,p4,p1D,p2D,p1L,p3L
#             p = l.nbrs[site];
#             sp = sb[p];
#             # some common factors in the matrix elements
#             sPw12 = simplePower(sp[1] + sp[2]);
#             sPw34 = simplePower(sp[3] + sp[4]);
#             sPw13 = simplePower(sp[1] + sp[3]);
#             sPw24 = simplePower(sp[2] + sp[4]);
#             sPw14 = simplePower(sp[1] + sp[4]);
#             sPw23 = simplePower(sp[2] + sp[3]);
#             sPw1234 = simplePower(sp[1] + sp[2] + sp[3] + sp[4]);
#             sPw1D2D = simplePower(sp[5] + sp[6]);
#             sPw1L3L = simplePower(sp[7] + sp[8]);
#
#             # some conditions
#             c12 = ((sPw12 == -1) && (Knz || J1nz));
#             c13 = ((sPw13 == -1) && (Knz || J1nz));
#             c14 = ((sPw14 == -1) && (Knz || J2nz));
#             c23 = ((sPw23 == -1) && (Knz || J2nz));
#             c1234 = ((sp[1] + sp[2] + sp[3] + sp[4] == 2) && Knz);
#
#
#             # contribute to the diagonal matrix element
#             # -----------------------------------------
#             if J1nz
#                 Hbb += 0.25*J1*( sPw12 + sPw13 );
#             end;
#             if J2nz
#                 Hbb += 0.25*J2*( sPw14 + sPw23 );
#             end;
#             if Knz
#                 Hbb += 0.125*K*sPw1234;
#             end;
#
#             # compute off diagonal matrix elements
#             # ------------------------------------
#
#             # the 12 term
#             if c12
#                 # the bra
#                 a = XiXj(b,p[1],p[2]);
#                 # the rep and transaltion of the bra
#                 aRep,lx,ly = representative(a,l);
#                 # search for this rep in the basis
#                 aRepIndex = basisIndex(aRep,basis);
#                 if aRepIndex != 0
#                     # the matrix element
#                     Hsite = (0.5*J1+0.125*K*(sPw34-sPw1234+2*sPw1D2D))*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);
#
#                     push!(I,aRepIndex);
#                     push!(M,Hsite);
#                     nnz += 1;
#                 end;
#             end;
#
#             # the 13 term
#             if c13
#                 # the bra
#                 a = XiXj(b,p[1],p[3]);
#                 # the rep and transaltion of the bra
#                 aRep,lx,ly = representative(a,l);
#                 # search for this rep in the basis
#                 aRepIndex = basisIndex(aRep,basis);
#                 if aRepIndex != 0
#                     # the matrix element
#                     Hsite = (0.5*J1+0.125*K*(sPw24-sPw1234+2*sPw1L3L))*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);
#
#                     push!(I,aRepIndex);
#                     push!(M,Hsite);
#                     nnz += 1;
#                 end;
#             end;
#
#             # the 14 term
#             if c14
#                 # the bra
#                 a = XiXj(b,p[1],p[4]);
#                 # the rep and transaltion of the bra
#                 aRep,lx,ly = representative(a,l);
#                 # search for this rep in the basis
#                 aRepIndex = basisIndex(aRep,basis);
#                 if aRepIndex != 0
#                     # the matrix element
#                     Hsite = (0.5*J2-0.25*K*sPw23)*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);
#
#                     push!(I,aRepIndex);
#                     push!(M,Hsite);
#                     nnz += 1;
#                 end;
#             end;
#
#             # the 23 term
#             if c23
#                 # the bra
#                 a = XiXj(b,p[2],p[3]);
#                 # the rep and transaltion of the bra
#                 aRep,lx,ly = representative(a,l);
#                 # search for this rep in the basis
#                 aRepIndex = basisIndex(aRep,basis);
#                 if aRepIndex != 0
#                     # the matrix element
#                     Hsite= (0.5*J2-0.25*K*sPw14)*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);
#
#                     push!(I,aRepIndex);
#                     push!(M,Hsite);
#                     nnz += 1;
#                 end;
#             end;
#
#             # the 1234 term
#             if c1234
#                 # the bra
#                 a = XiXjXkXl(b,p[1],p[2],p[3],p[4]);
#                 # the rep and transaltion of the bra
#                 aRep,lx,ly = representative(a,l);
#                 # search for this rep in the basis
#                 aRepIndex = basisIndex(aRep,basis);
#                 if aRepIndex != 0
#                     # the matrix element
#                     Hsite = 0.125*K*(2-sPw12-sPw34-sPw13-sPw24+sPw14+sPw23)*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);
#
#                     push!(I,aRepIndex);
#                     push!(M,Hsite);
#                     nnz += 1;
#                 end;
#             end;
#
#         end;
#
#         # push diagonal matrix element to list of matrix elements
#         push!(I,bIndex);
#         push!(M,Hbb);
#         nnz += 1;
#
#         # CSC formatting (turns out this is unnecessary)
#         # sortTwo!(I,M,Int64(Jpointers[bIndex]),Int64(Jpointers[bIndex])+nnz-1);
#
#
#     end;
#
#     # CSC formatting
#     Jpointers[end] = Int32(length(I)+1);
#
#     H::SparseMatrixCSC{Complex128,Int32} = SparseMatrixCSC{Complex128,Int32}(basis.dim, basis.dim, Jpointers, I, M);
#     return H;
# end;

# function for building the dense Hamiltonian
function constructDenseHam(basis::SzkxkyBasis,c::couplings,s::sector,l::lattice)

    # lattice
    Lx::Int64 = l.Lx;
    Ly::Int64 = l.Ly;
    N::Int64 = l.N;
    # couplings
    J1::Float64 = c.J1;
    J2::Float64 = c.J2;
    K::Float64 = c.K;
    # momentum
    kx::Float64 = s.kx;
    ky::Float64 = s.ky;

    # non-zero couplings
    J1nz::Bool = (J1 != 0.0);
    J2nz::Bool = (J2 != 0.0);
    Knz::Bool = (K != 0.0);

    # basis type
    bType::Type = eltype(basis);

    # allocate memory once before the loops
    # -------------------------------------
    H::Matrix{Complex128} = zeros(Complex128,basis.dim,basis.dim);
    b::bType = 0;
    sb::Array{bType,1} = Array{bType,1}(l.N);
    nb::Float64 = 0.0;
    Hbb::Complex128 = 0.0+0.0im;
    p::Array{Int64,1} = Array{Int64,1}(length(l.nbrs[1]));
    sp::Array{bType,1} = Array{UInt64,1}(length(p));
    sPw12::Int64 = sPw34::Int64 = sPw13::Int64 = sPw24::Int64 = sPw14::Int64 = sPw23::Int64 = sPw1234::Int64 = sPw1D2D::Int64 = sPw1L3L::Int64 = 0;
    c12::Bool = c13::Bool = c14::Bool = c23::Bool = c1234::Bool = true;
    a::bType = convert(bType,0);
    aRep::bType = convert(bType,0);
    lx::Int64 = ly::Int64 = 0;
    aRepIndex::Int32 = Int32(0);
    Hsite::Complex128 = 0.0+0.0im;
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

        # initialize the diagonal matrix element
        Hbb = 0.0 + 0.0im;

        # loop over lattice sites
        for site::Int64 in 1:N

            # spins at nearby sites p1,p2,p3,p4,p1D,p2D,p1L,p3L
            p = l.nbrs[site];
            sp = sb[p];
            # some common factors in the matrix elements
            sPw12 = simplePower(sp[1] + sp[2]);
            sPw34 = simplePower(sp[3] + sp[4]);
            sPw13 = simplePower(sp[1] + sp[3]);
            sPw24 = simplePower(sp[2] + sp[4]);
            sPw14 = simplePower(sp[1] + sp[4]);
            sPw23 = simplePower(sp[2] + sp[3]);
            sPw1234 = simplePower(sp[1] + sp[2] + sp[3] + sp[4]);
            sPw1D2D = simplePower(sp[5] + sp[6]);
            sPw1L3L = simplePower(sp[7] + sp[8]);

            # some conditions
            c12 = ((sPw12 == -1) && (Knz || J1nz));
            c13 = ((sPw13 == -1) && (Knz || J1nz));
            c14 = ((sPw14 == -1) && (Knz || J2nz));
            c23 = ((sPw23 == -1) && (Knz || J2nz));
            c1234 = ((sp[1] + sp[2] + sp[3] + sp[4] == 2) && Knz);


            # contribute to the diagonal matrix element
            # -----------------------------------------
            if J1nz
                Hbb += 0.25*J1*( sPw12 + sPw13 );
            end;
            if J2nz
                Hbb += 0.25*J2*( sPw14 + sPw23 );
            end;
            if Knz
                Hbb += 0.125*K*sPw1234;
            end;

            # compute off diagonal matrix elements
            # ------------------------------------

            # the 12 term
            if c12
                # the bra
                a = XiXj(b,p[1],p[2]);
                # the rep and transaltion of the bra
                aRep,lx,ly = representative(a,l);
                # search for this rep in the basis
                aRepIndex = basisIndex(aRep,basis);
                if aRepIndex != 0
                    # the matrix element
                    H[aRepIndex,bIndex] +=  (0.5*J1+0.125*K*(sPw34-sPw1234+2*sPw1D2D))*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);

                end;
            end;

            # the 13 term
            if c13
                # the bra
                a = XiXj(b,p[1],p[3]);
                # the rep and transaltion of the bra
                aRep,lx,ly = representative(a,l);
                # search for this rep in the basis
                aRepIndex = basisIndex(aRep,basis);
                if aRepIndex != 0
                    # the matrix element
                    H[aRepIndex,bIndex] += (0.5*J1+0.125*K*(sPw24-sPw1234+2*sPw1L3L))*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);

                end;
            end;

            # the 14 term
            if c14
                # the bra
                a = XiXj(b,p[1],p[4]);
                # the rep and transaltion of the bra
                aRep,lx,ly = representative(a,l);
                # search for this rep in the basis
                aRepIndex = basisIndex(aRep,basis);
                if aRepIndex != 0
                    # the matrix element
                   H[aRepIndex,bIndex] += (0.5*J2-0.25*K*sPw23)*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);

                end;
            end;

            # the 23 term
            if c23
                # the bra
                a = XiXj(b,p[2],p[3]);
                # the rep and transaltion of the bra
                aRep,lx,ly = representative(a,l);
                # search for this rep in the basis
                aRepIndex = basisIndex(aRep,basis);
                if aRepIndex != 0
                    # the matrix element
                    H[aRepIndex,bIndex] += (0.5*J2-0.25*K*sPw14)*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);

                end;
            end;

            # the 1234 term
            if c1234
                # the bra
                a = XiXjXkXl(b,p[1],p[2],p[3],p[4]);
                # the rep and transaltion of the bra
                aRep,lx,ly = representative(a,l);
                # search for this rep in the basis
                aRepIndex = basisIndex(aRep,basis);
                if aRepIndex != 0
                    # the matrix element
                    H[aRepIndex,bIndex] += 0.125*K*(2-sPw12-sPw34-sPw13-sPw24+sPw14+sPw23)*exp(-1.0im*(kx*lx+ky*ly))*sqrt(basis.n[aRepIndex]/nb);

                end;
            end;

        end;

        # push diagonal matrix element to list of matrix elements
        H[bIndex,bIndex] += Hbb;

    end;

    return H;
end;
