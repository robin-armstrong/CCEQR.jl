using LinearAlgebra

### Fills the (:, p:q) block of V and the (1:q, p:q) block of T for compact WY form using the Schreiber-van Loan algorithm.

function update_q!(T::AbstractMatrix{R}, V::Matrix{R}, F::Matrix{R}, tau::Vector{R}, p::Int, q::Int) where R <: Real
    # filling the (:, p:q) block of V
    
    V[p:end, p:q] = F[:, 1:(q-p+1)]

    for i = p:q
        V[i, i] = 1.
        fill!(view(V, 1:(i-1), i), 0.)
    end
    
    # filling the (p:q, p:q) block of T
    
    for s = 1:(q-p+1)
        i       = p+s-1
        T[i, i] = tau[s]
        
        if i > p
            c      = view(T, p:(i-1), i)
            V_prev = view(V, :, p:(i-1))
            V_new  = view(V, :, i)
            T_prev = UpperTriangular(view(T, p:(i-1), p:(i-1)))

            mul!(c, V_prev', V_new, -tau[s], 0.)
            lmul!(T_prev, c)
        end
    end
    
    # filling the (1:(p-1), p:q) block of T
    
    if p > 1
        T11 = UpperTriangular(view(T, 1:(p-1), 1:(p-1)))
        T22 = UpperTriangular(view(T, p:q, p:q))
        T12 = view(T, 1:(p-1), p:q)
        V1  = view(V, :, 1:(p-1))
        V2  = view(V, :, p:q)

        mul!(T12, V1', V2, -1., 0.)
        lmul!(T11, T12)
        rmul!(T12, T22)
    end
end

### Applies Qt to the (p:q, r:s) block of A, where Qt is defined by the p-th through q-th Householder reflectors
### in the compact WY form given by V and T.

function apply_qt!(A::Matrix{R}, V::Matrix{R}, T::Matrix{R}, p::Int, q::Int, r::Int, s::Int) where R <: Real
    m, n  = size(A)
    A1    = view(A, p:q, r:s)
    A2    = view(A, (q+1):m, r:s)
    T_sub = view(T, p:q, p:q)
    V1    = view(V, p:q, p:q)
    V2    = view(V, (q+1):m, p:q)
    W     = Matrix{R}(A[p:q, r:s]')

    rmul!(W, LowerTriangular(V1))
    mul!(W, A2', V2, 1., 1.)
    rmul!(W, UpperTriangular(T_sub))
    mul!(A1, V1, W', -1., 1.)
    mul!(A2, V2, W', -1., 1.)
end

### Loops through the (:, j_start:j_end) block of A and moves all columns to the front that have squared norm exceeding threshold
### thresh. Modifies jpvt and gamma accordingly. Returns the number of columns that passed the threshold, and the maximum norm
### among columns that did not pass the threshold.

function threshold_reblock!(A::Matrix{R}, jpvt::Vector{T}, j_start::Int, j_end::Int, gamma::Vector{R}, thresh::Real) where {R <: Real, T <: Int}
    m, n   = size(A)
    tmpcol = zeros(m)
    blk    = 0
    maxut  = 0.

    for j = j_start:j_end
        if gamma[j] > thresh
            blk += 1
            p    = j_start+blk-1
            
            tmp     = jpvt[p]
            jpvt[p] = jpvt[j]
            jpvt[j] = tmp

            tmp      = gamma[p]
            gamma[p] = gamma[j]
            gamma[j] = tmp

            tmpcol[:] = A[:,p]
            A[:,p]    = A[:,j]
            A[:,j]    = tmpcol
        else
            maxut = max(maxut, gamma[j])
        end
    end

    return blk, maxut
end

### Loops through the (:, j_start:j_end) block of A and moves r columns to the front, corresponding to the columns
### with the highest value of gamma. Returns the (r+1)th largest gamma value in the block, i.e., the largest gamma
### value *not* included in this set. Modifies jpvt and gamma accordingly.

function order_reblock!(A::Matrix{R}, jpvt::Vector{T}, j_start::Int, j_end::Int, gamma::Vector{R}, r::Int) where {R <: Real, T <: Int}
    m, n   = size(A)
    tmpcol = zeros(m)
    gview  = view(gamma, j_start:j_end)

    if r < j_end - j_start
        samp = partialsortperm(gview, 1:(r+1), rev = true)
        g    = gview[samp[r+1]]
    else
        samp = sortperm(gview, rev = true)
        g    = 0.
    end

    sort!(view(samp, 1:r))

    for i = 1:r
        p = j_start+samp[i]-1
        j = j_start+i-1

        tmp     = jpvt[p]
        jpvt[p] = jpvt[j]
        jpvt[j] = tmp

        tmp      = gamma[p]
        gamma[p] = gamma[j]
        gamma[j] = tmp

        tmpcol[:] = A[:,p]
        A[:,p]    = A[:,j]
        A[:,j]    = tmpcol
    end

    return g
end

### Applies column permutation "perm" to the (:, j0:(j0+length(perm)-1)) block of A, modifying jpvt and gamma accordingly.

function swap_cols!(A::Matrix{R}, jpvt::Vector{T}, gamma::Vector{R}, j0::Int, perm::Vector{T}) where {R <: Real, T <: Int}
    n = size(A, 2)
    b = length(perm)

    permute!(view(jpvt, j0:(j0+b-1)), perm)
    permute!(view(gamma, j0:(j0+b-1)), perm)
    Base.permutecols!!(view(A, :, j0:(j0+b-1)), perm)
end
