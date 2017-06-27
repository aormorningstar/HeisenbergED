# sparseHermitian.jl
# sparse Hermitian matrix type and matrix vector multiplication
# Alan Morningstar
# June 2017

immutable sparseHermitian

    # NOTE: diagonal elements here should have been divided by 2, only the upper triangle of the off diagonal elements are included

    # number of rows in the matrix
    nRows::Int64;
    # number of columns in the matrix
    nCols::Int64;
    # indices of the first non-zero value of each column
    colPntrs::Array{Int32,1};
    # row indices of non-zero elements
    rowIndcs::Array{Int32,1};
    # non-zero matrix element values
    nzVals::Array{Complex128,1};

end;

function matVecSparseHermitian!(M::sparseHermitian,y::AbstractVector,x::AbstractVector)

    # NOTE: this assumes the sparse Hermitian matrix has been built properly, with diagonal matrix elements divided by 2, and only the upper triangle of the off-diagonal elements included

    # clear output vector
    y .= 0.0 + 0.0im;

    # fill output vector
    for j::Int64 in 1:M.nCols
        for i::Int64 in M.colPntrs[j]:(M.colPntrs[j+1]-1)
            y[M.rowIndcs[i]] += x[j]*M.nzVals[i];
            y[j] += x[M.rowIndcs[i]]*conj(M.nzVals[i]);
        end;
    end;

    return y::AbstractVector;

end;
