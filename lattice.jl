# lattice.jl
# defining the lattice connectivity for rectangular 2D lattice with periodic boundaries
# Alan Morningstar
# May 2017


# convert single integer site index to xy indexing
function xyIndex(site::Int64,Lx::Int64,Ly::Int64)
    # site - number of the site for which x,y indices are computed

    y::Int,x::Int = divrem(site-1,Lx);

    return x,y;
end;


# convert xy site indexing to single integer index
function siteIndex(xy::Tuple{Int,Int},Lx::Int64,Ly::Int64)
    # xy - tuple containing x,y indexing of a site

    x = xy[1];
    y = xy[2];

    # move xy into primary lattice cell
    while x > (Lx-1)
        x -= Lx;
    end;
    while x < 0
        x += Lx;
    end;
    while y > (Ly-1)
        y -= Ly;
    end;
    while y < 0
        y += Ly;
    end;

    # map to integer index
    site::Int64 = y*Lx + x + 1;

    return site;
end;


# find indices of neighbors specified by a translation vectors T
function neighbors(site::Int64,T::Array{Tuple{Int64,Int64},1},Lx::Int64,Ly::Int64)
    # use xy indexing just to compute neighbors
    x::Int,y::Int = xyIndex(site,Lx,Ly);

    # neighbor in xy indexing
    neighborxy::Array{Tuple{Int,Int},1} = [(x+t[1],y+t[2]) for t in T];

    # back to single integer indexing
    neighborIndex::Array{Int,1} = [siteIndex(nxy,Lx,Ly) for nxy in neighborxy];

    # return neighbor index
    return neighborIndex;
end;


#-- define the lattice

# plaquette sites corresponding to a given site (p1) in the bottom left of the plaquette
# p3L indicates p3 on the plaquette to the left, similarly for p1D on the plaquette in the downwards direction
#
# p3L--p3--p4
#   |  |Plq|
# p1L--p1--p2
#  |  |   |
#   p1D p2D


# container for lattice properties
immutable lattice
    # size of lattice, number of sites
    Lx::Int64;
    Ly::Int64;
    N::Int64;
    nbrs::Array{Array{Int64,1},1};
    TxMask1::UInt64;
    TxMask2::UInt64;
    TyMask1::UInt64;
    TyMask2::UInt64;
    LxMinus1::Int64;
    NMinusLx::Int64;

    function lattice(Lx::Int64,Ly::Int64,neighborVectors::Array{Tuple{Int64,Int64},1})
        neighborsList = [neighbors(site,neighborVectors,Lx,Ly) for site in 1:Lx*Ly];

        # rough code, needs cleanup
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

        new(Lx,Ly,Lx*Ly,neighborsList,TxMask1,TxMask2,TyMask1,TyMask2,Lx-1,Lx*(Ly-1));
    end;
end;


# x translation operator
function Tx(b::UInt64,l::lattice)
    return ((b & l.TxMask1) >> (l.LxMinus1) | ((b & l.TxMask2) << 1));
end;


# y translation operator
function Ty(b::UInt64,l::lattice)
    return ((b & l.TyMask1) >> (l.NMinusLx) | ((b & l.TyMask2) << l.Lx));
end;


# find related representative state and what translation relates the two states
function representative(b::UInt64,l::lattice)
    # rep. state
    rep::UInt64 = b;
    # translated state
    Tb::UInt64 = b;
    # translations to get to rep. state
    lx::Int64 = 0;
    ly::Int64 = 0;

    # perform all translations
    for y::Int in 0:l.Ly-1
        for x::Int in 0:l.Lx-1

            if Tb < rep
                # then this is a better representative
                rep = Tb;
                lx = x;
                ly = y;
            end;

            Tb = Tx(Tb,l);
        end;
        Tb = Ty(Tb,l);
    end;

    return rep::UInt64,lx::Int64,ly::Int64;
end;
