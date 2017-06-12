# utils.jl
# useful functions for ED project
# Alan Morningstar
# May 2017


# set a bit to 1
function setBit(i::UInt64,bit::Int64)
    return (i | (1 << UInt64(bit-1)))::UInt64;
end;


# set a bit to 0
function clearBit(i::UInt64,bit::Int64)
    return (i & (~(UInt64(1) << UInt64(bit-1))))::UInt64;
end;


# toggle a bit
function toggleBit(i::UInt64,bit::Int64)
    return (i $ (UInt64(1) << UInt64(bit-1)))::UInt64;
end;


# read a bit
function readBit(i::UInt64,bit::Int64)
    return ((i >> UInt64(bit-1)) & UInt64(1))::UInt64;
end;


# swap two bits
function swapBits(i::UInt64,bit1::Int64,bit2::Int64)
    if readBit(i,bit1) == readBit(i,bit2)
        return i::UInt64;
    else
        return (i $ ( (UInt64(1)<<UInt64(bit1-1)) | (UInt64(1)<<UInt64(bit2-1)) ))::UInt64;
    end;
end;


# count number of unique elements in an array
function numUnique(Tbs::Array{UInt64,1})
    # num unique elements
    num::Int64 = 0;
    # scan over elements looking in previous elements
    for i in 1:length(Tbs)
        if !(Tbs[i] in Tbs[1:i-1])
            num += 1;
        end;
    end;

    return num::Int64;
end;


# function for flipping two spins of a basis state
function XiXj(b::UInt64,i::Int64,j::Int64)
    # b - integer rep of state
    # i,j - sites of spin flips
    return b $ ( (UInt64(1)<<UInt64(i-1)) | (UInt64(1)<<UInt64(j-1)) );
end;


# function for flipping four spins of a basis state
function XiXjXkXl(b::UInt64,i::Int,j::Int,k::Int,l::Int)
    # b - integer rep of state
    # i,j,k,l - sites of spin flips
    return b $ ( (UInt64(1)<<UInt64(i-1)) | (UInt64(1)<<UInt64(j-1)) | (UInt64(1)<<UInt64(k-1)) | (UInt64(1)<<UInt64(l-1)) );
end;


# function for computing (-1)^(positive integer)
function simplePower(i::UInt64)
    if readBit(i,1) == UInt64(0)
        return 1::Int64;
    else
        return -1::Int64;
    end;
end;

# function for sorting two arrays in place based on the first of the arrays
function sortTwo!(I::Array{Int32,1},M::Array{Complex128,1},low::Int64,high::Int64)
    if low < high
        p::Int64 = partition!(I,M,low,high);
        sortTwo!(I,M,low,p);
        sortTwo!(I,M,p+1,high);
    end;
end;
function partition!(I::Array{Int32,1},M::Array{Complex128,1},low::Int64,high::Int64)
    piv::Int32 = I[low];
    i::Int64 = low;
    j::Int64 = high;

    while true
        while I[i] < piv
            i += 1;
        end;
        while I[j] > piv
            j -= 1;
        end;
        if i >= j
          return j;
        end;
        I[i],I[j],M[i],M[j] = I[j],I[i],M[j],M[i];
        i += 1;
        j -= 1;
    end;
end;
