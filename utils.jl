# utils.jl
# useful functions for ED project
# Alan Morningstar
# May 2017


# set a bit to 1
function setBit{I<:Integer}(i::I,bit::Int64)
    return (i | (1 << (bit-1)))::I;
end;


# set a bit to 0
function clearBit{I<:Integer}(i::I,bit::Int64)
    return convert(I,(Int64(i) & (~(1 << (bit-1)))));
end;


# toggle a bit
function toggleBit{I<:Integer}(i::I,bit::Int64)
    return xor(i , (1 << (bit-1)))::I;
end;


# read a bit
function readBit{I<:Integer}(i::I,bit::Int64)
    return Int64((i >> (bit-1)) & 1);
end;


# swap two bits
function swapBits{I<:Integer}(i::I,bit1::Int64,bit2::Int64)
    if readBit(i,bit1) == readBit(i,bit2)
        return i::I;
    else
        return xor(i , ( (1<<(bit1-1)) | (1<<(bit2-1)) ))::I;
    end;
end;


# count number of unique elements in an array of Ints
function numUnique!{I<:Real}(Tbs::Array{I,1})
    # first sort the list in place
    sort!(Tbs);
    # num unique elements
    num::Int64 = 1;
    # scan over elements looking in previous elements
    for i::Int64 in 2:length(Tbs)
        if Tbs[i] > Tbs[i-1]
            num += 1;
        end;
    end;

    return num::Int64;
end;


# function for flipping two spins of a basis state
function XiXj{I<:Integer}(b::I,i::Int64,j::Int64)
    # b - integer rep of state
    # i,j - sites of spin flips
    return xor(b , ( (1<<(i-1)) | (1<<(j-1)) ))::I;
end;


# function for flipping four spins of a basis state
function XiXjXkXl{I<:Integer}(b::I,i::Int64,j::Int64,k::Int64,l::Int64)
    # b - integer rep of state
    # i,j,k,l - sites of spin flips
    return xor(b , ( (1<<(i-1)) | (1<<(j-1)) | (1<<(k-1)) | (1<<(l-1)) ))::I;
end;


# function for computing (-1)^(positive integer)
function simplePower{I<:Integer}(i::I)
    if readBit(i,1) == 0
        return 1::Int64;
    else
        return -1::Int64;
    end;
end;


# function for sorting multiple lists in place based on the values of another list which is also sorted
function sortTwo!{TI<:Real}(low::Int64,high::Int64,I::Vector{TI},M...)
    if low < high
        p::Int64 = partition!(low,high,I,M...);
        sortTwo!(low,p,I,M...);
        sortTwo!(p+1,high,I,M...);
    end;
end;
function partition!{TI<:Real}(low::Int64,high::Int64,I::Vector{TI},M...)
    piv::TI = I[low];
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
        I[i],I[j] = I[j],I[i];
        for m in M
            m[i],m[j] = m[j],m[i]
        end;
        i += 1;
        j -= 1;
    end;
end;


# find the index of the first instance of an element in a set (vector), if it's not in the set then append the element and return it's index
function appendSet!{T}(elem::T,set::Vector{T})

    for (index::Int64,item::T) in enumerate(set)
        if elem == item
            return index::Int64;
        end;
    end;

    # if it wasn't found, append the element
    push!(set,elem);

    return length(set)::Int64;

end;
