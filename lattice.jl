# lattice.jl
# defining the lattice connectivity for rectangular 2D lattice with
# periodic boundaries
# Alan Morningstar
# May 2017


# convert single integer site index to xy indexing
# NOTE: lattice sites are labelled like the following example
#
# integer labelling
# -----------------
# 34
# 12
#
# xy labelling
# ------------
# (0,1)(1,1)
# (0,0)(1,0)
#
function xyIndex(site::Int64,Lx::Int64,Ly::Int64)

    y::Int64,x::Int64 = divrem(site-1,Lx)

    return x,y

end


# convert xy site indexing to single integer index
function siteIndex(xy::Tuple{Int64,Int64},Lx::Int64,Ly::Int64)

    x::Int64 = xy[1]
    y::Int64 = xy[2]

    # move xy into primary lattice cell
    while x > (Lx-1)
        x -= Lx
    end
    while x < 0
        x += Lx
    end
    while y > (Ly-1)
        y -= Ly
    end
    while y < 0
        y += Ly
    end

    # map to integer index
    site::Int64 = y*Lx + x + 1

    return site

end


# find indices of neighbors specified by a translation vectors T.
# takes in tuple of vectors pointing to so-called neighbors (this is loosely
# defined and can be used to collect---for each site---a list of 'neighboring'
# sites that are connected to it by the Hamiltonian)
# ex: if the Hamiltonian couples plaquette spins and next nearest neighbors,
#     then include those in this list so they are easily accessed by other
#     customizable functions (like building the Hamiltonian)
function neighbors(site::Int64,T::Array{Tuple{Int64,Int64},1},Lx::Int64,Ly::Int64)

    # use xy indexing just to compute neighbors
    x::Int64,y::Int64 = xyIndex(site,Lx,Ly)

    # neighbor in xy indexing
    neighborxy::Array{Tuple{Int64,Int64},1} = [(x+t[1],y+t[2]) for t in T]

    # back to single integer indexing
    neighborIndex::Array{Int64,1} = [siteIndex(nxy,Lx,Ly) for nxy in neighborxy]

    # return neighbor index
    return neighborIndex

end


# NOTE: plaquette sites corresponding to a given site (p1) in the bottom left
# of the plaquette.
# p3L indicates p3 on the plaquette to the left, similarly for p1D on the
# plaquette in the downwards direction
#
# p3L--p3--p4
#   |  |Plq|
# p1L--p1--p2
#  |  |   |
#   p1D p2D
#
# this is useful for J1-J2-K model


# container for lattice properties
immutable lattice

    # size of lattice, number of sites
    Lx::Int64
    Ly::Int64
    N::Int64

    # neighbors of each site detailed in diagram above
    nbrs::Array{Array{Int64,1},1}

    # useful constants for translation and inversion operators
    TxMask1::UInt64
    TxMask2::UInt64
    TyMask1::UInt64
    TyMask2::UInt64
    LxMinus1::Int64
    NMinusLx::Int64
    ZMask::UInt64

    function lattice(Lx::Int64,Ly::Int64,neighborVectors::Array{Tuple{Int64,Int64},1})

        neighborsList = [neighbors(site,neighborVectors,Lx,Ly) for site in 1:Lx*Ly]

        # rough code, needs cleanup.
        # just makes the useful masks
        # ------------------------------------------
        TxMx1::UInt64 = UInt64(1)<<(Lx-1)
        TxMx2::UInt64 = TxMx1 - UInt64(1)
        TxMask1::UInt64,TxMask2::UInt64 = TxMx1,TxMx2

        for y in 1:Ly
            TxMask1 = TxMask1 | (TxMx1 << y*Lx)
            TxMask2 = TxMask2 | (TxMx2 << y*Lx)
        end

        TyMask1 = ((UInt64(1)<<Lx)-1) << (Lx * (Ly-1))
        TyMask2 = (UInt64(1)<< (Lx * (Ly-1))) - 1
        # -------------------------------------------

        new(Lx,Ly,Lx*Ly,neighborsList,TxMask1,TxMask2,TyMask1,TyMask2,Lx-1,Lx*(Ly-1),UInt64(2^(Lx*Ly)-1))

    end

end


# x translation operator, shifts spin values to the right on the lattice
function Tx{I<:Integer}(b::I,l::lattice)
    return ((b & l.TxMask1) >> (l.LxMinus1) | ((b & l.TxMask2) << 1))
end


# y translation operator, shifts spin values up on the lattice
function Ty{I<:Integer}(b::I,l::lattice)
    return ((b & l.TyMask1) >> (l.NMinusLx) | ((b & l.TyMask2) << l.Lx))
end


# spin flip operator, flips all spins (bits)
function Z{I<:Integer}(b::I,l::lattice)
    return ( l.ZMask $ b )
end


# find the related representative state and what translation and spin flip
# relate the two states
function representative{I<:Integer}(b::I,l::lattice)

    # rep. state
    rep::I = b
    # transformed states
    Tb::I = b
    ZTb::I = Z(b,l)
    # translations to get to rep. state
    lx::Int64 = 0
    ly::Int64 = 0
    # inversion to get to rep. state
    g::Int64 = 0;

    # perform all transformations
    for y::Int64 in 0:l.Ly-1
        for x::Int64 in 0:l.Lx-1

            if Tb < rep
                # then this is a better representative
                rep = Tb
                lx = x
                ly = y
                g = 0
            end

            # now invert the spins
            ZTb = Z(Tb,l)

            if ZTb < rep
                # then this is a better representative
                rep = ZTb
                lx = x
                ly = y
                g = 1
            end


            Tb = Tx(Tb,l)

        end

        Tb = Ty(Tb,l)

    end

    # return lx, ly, and g are the powers of their respective operators
    # Tx,Ty,Z which---when applied to the state b---yield the representative
    return rep::I,lx::Int64,ly::Int64,g::Int64

end
