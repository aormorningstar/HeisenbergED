# sparseHermitian.jl
# sparse Hermitian matrix type and matrix vector multiplication
# Alan Morningstar
# June 2017

struct sparseHermitian{TypeI,TypeM}

    # NOTE: diagonal elements here should have been divided by 2, only the upper triangle of the off diagonal elements are included

    # number of rows in the matrix
    nRows::Int64;
    # number of columns in the matrix
    nCols::Int64;
    # indices of the first non-zero value of each column
    colPntrs::Vector{TypeI};
    # row indices of non-zero elements
    rowIndcs::Vector{TypeI};
    # non-zero matrix element values
    nzVals::Vector{TypeM};

end;

function Base.A_mul_B!(y::AbstractVector,M::sparseHermitian,x::AbstractVector)

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

Base.size(M::sparseHermitian) = (M.nRows, M.nCols);

Base.length(M::sparseHermitian) = reduce(*, size(M))

Base.ndims(M::sparseHermitian) = 2;

Base.eltype(M::sparseHermitian) = eltype(M.nzVals);

Base.issymmetric(M::sparseHermitian) = false;

Base.ishermitian(M::sparseHermitian) = true;

Base.nnz(M::sparseHermitian) = length(M.nzVals);

Base.:*(M::sparseHermitian,x::AbstractVector) = Base.A_mul_B!(similar(x,promote_type(Base.eltype(M),eltype(x)),size(M,1)),M,x);

Base.Ac_mul_B!(y::AbstractVector,M::sparseHermitian,x::AbstractVector) = Base.A_mul_B!(y,M,x);
