T = Union{SubArray{ComplexF64}, Vector{ComplexF64}}

"""
Rearrange 1d-arrays's values in bit reverse order in-place
"""
function _bit_reverse_permutation_1d!(arr::T)
    n = size(arr, 1)
    j = 1
    for i in 2:n
        bit = n >> 1
        while j > bit
            j -= bit
            bit >>= 1
        end
        j += bit
        if i < j
            arr[i], arr[j] = arr[j], arr[i]
        end
    end
end

function _fft1_pow2_process_block2!(arr::T, block::Int)
    arr[block + 1] += arr[block + 2]
    arr[block + 2] = arr[block + 1] - 2*arr[block + 2]
end

function _fft1_pow2_process_block!(arr::T, half::Int, block::Int, roots::Vector{ComplexF64})
    @fastmath @inbounds @simd for i in 1:half
        arr[block + i + half] *= roots[i]
        arr[block + i] += arr[block + i + half]
        arr[block + i + half] = arr[block + i] - 2*arr[block + i + half]
    end
end

function _fft1_pow2!(arr::T, roots::Union{Vector{ComplexF64}, Nothing}=nothing)
    typeof(arr) <: AbstractArray || error("got non-array")
    ndims(arr) == 1              || error("got more than 1-d array")
    ispow2(size(arr, 1))         || error("array's length is not power of 2")

    set_zero_subnormals(true)
    n = size(arr, 1)
    _bit_reverse_permutation_1d!(arr)
    if roots === nothing
        roots = Array{ComplexF64}(undef, n)
    end
    if n >= 2
        @fastmath @inbounds _fft1_pow2_process_block2!.(Ref(arr), 0:2:n-1)
    end
    len = 4
    while len <= n
        half = len >> 1
        @fastmath @inbounds @simd for i in 0:half - 1
            roots[i + 1] = cispi(-2 * i / len)
        end
        @fastmath @inbounds _fft1_pow2_process_block!.(Ref(arr), half, 0:len:n-1, Ref(roots))
        len <<= 1
    end
    set_zero_subnormals(false)
end

function _fft1_pow2(arr::Array, axis::Int)
    arr = ComplexF64.(arr)  # make a complex copy
    dims = ndims(arr)
    len = size(arr, axis)
    if dims == 1
        _fft1_pow2!(arr)
    else
        idxs = tuple(deleteat!(collect(1:dims), axis)...)
        arr_perm = axis != 1 ? permutedims(arr, (axis, idxs...)) : arr
        slices = eachslice(arr_perm, dims=(2:dims...,))
        roots = Array{ComplexF64}(undef, len)
        _fft1_pow2!.(slices, Ref(roots))
        if axis != 1
            permutedims!(arr, arr_perm, (2:axis..., 1, axis+1:dims...))
        end
    end
    return arr
end

"""
Rearrange 2d-arrays's columns and rows in bit reverse order in-place
"""
function _bit_reverse_permutation_2d!(arr::Matrix{ComplexF64}, axis::Int=1)
    m, n = size(arr)
    j = 1
    for i in 2:(axis==1 ? n : m)
        bit = n >> 1
        while j > bit
            j -= bit
            bit >>= 1
        end
        j += bit
        if i < j
            for k in 1:(axis==1 ? m : n)
                if axis == 1
                    arr[k, i], arr[k, j] = arr[k, j], arr[k, i]
                else
                    arr[i, k], arr[j, k] = arr[j, k], arr[i, k]
                end
            end
        end
    end
end

function _fft2_pow2_process_block2!(arr::Matrix{ComplexF64}, block1::Int, block2::Int)
    @fastmath @inbounds begin
        u11 = arr[block1 + 1, block2 + 1]
        u21 = arr[block1 + 2, block2 + 1]
        u22 = arr[block1 + 2, block2 + 2]
        u12 = arr[block1 + 1, block2 + 2]

        arr[block1 + 1, block2 + 2] = u11 - u12 + u21 - u22  # 12
        arr[block1 + 1, block2 + 1] = u11 + u12 + u21 + u22  # 11
        arr[block1 + 2, block2 + 1] = u11 + u12 - u21 - u22  # 21
        arr[block1 + 2, block2 + 2] = u11 - u12 - u21 + u22  # 22
    end

    # arr[block1 + 1, block2 + 1] += arr[block1 + 2, block2 + 1] + arr[block1 + 2, block2 + 2] + arr[block1 + 1, block2 + 2]
    # u12 = arr[block1 + 1, block2 + 2]
    # u21 = arr[block1 + 2, block2 + 1]
    # arr[block1 + 1, block2 + 2] = arr[block1 + 1, block2 + 1] - 2 * (arr[block1 + 1, block2 + 2] + arr[block1 + 2, block2 + 2])
    # arr[block1 + 2, block2 + 1] = arr[block1 + 1, block2 + 1] - 2 * (arr[block1 + 2, block2 + 1] + arr[block1 + 2, block2 + 2])
    # arr[block1 + 2, block2 + 2] = arr[block1 + 1, block2 + 1] - 2 * (u12 + u21)
end

function _fft2_pow2_cooley_square!(arr::Matrix{ComplexF64})
    typeof(arr) <: AbstractArray || error("got non-array")
    ndims(arr) == 2              || error("got more than 2-d array")
    ispow2(size(arr, 1))         || ispow2(size(arr, 2)) || error("array's dimensions are not power of 2")

    set_zero_subnormals(true)

    n = size(arr, 1)

    _bit_reverse_permutation_2d!(arr)
    arr_perm = permutedims(arr, (2, 1))
    _bit_reverse_permutation_2d!(arr_perm)
    permutedims!(arr, arr_perm, (2, 1))

    if n >= 2
        @fastmath @inbounds @simd for block2 in 0:2:n-1
            @fastmath @inbounds @simd for block1 in 0:2:n-1
                @fastmath @inbounds _fft2_pow2_process_block2!(arr, block1, block2)
            end
        end
    end

    len = 4
    roots = Array{ComplexF64}(undef, n)
    while len <= n
        half = len >> 1
        @fastmath @inbounds @simd for i in 0:half - 1
            roots[i + 1] = cispi(-2 * i / len)
        end
        @fastmath @inbounds @simd for block2 in 0:len:n-1
            @fastmath @inbounds @simd for j in 1:half
                @fastmath @inbounds @simd for block1 in 0:len:n-1
                    @fastmath @inbounds @simd for i in 1:half
                        @fastmath @inbounds begin
                            u11 = arr[block1 + i, block2 + j]
                            u21 = arr[block1 + i + half, block2 + j] * roots[i]
                            u22 = arr[block1 + i + half, block2 + j + half] * roots[i] * roots[j]
                            u12 = arr[block1 + i, block2 + j + half] * roots[j]

                            arr[block1 + i, block2 + j] = u11 + u12 + u21 + u22
                            arr[block1 + i + half, block2 + j] = u11 + u12 - u21 - u22
                            arr[block1 + i + half, block2 + j + half] = u11 - u12 - u21 + u22
                            arr[block1 + i, block2 + j + half] = u11 - u12 + u21 - u22
                        end
                    end
                end
            end
        end
        len <<= 1
    end

    set_zero_subnormals(false)
end

function _fft2_pow2_cooley_rect!(arr::Matrix{ComplexF64})
    n, m = size(arr)
    if n == m
        _fft2_pow2_cooley_square!(arr)
        return
    end

    new_m = m >> 1
    tmp = fill(Matrix{ComplexF64}(undef, n, new_m), 2)

    tmp[1] = arr[:, 1:2:m]
    tmp[2] = arr[:, 2:2:m]

    _fft2_pow2_cooley_rect!.(tmp[1:2])

    wnp = cispi.(-2 * (0:new_m-1) / m)
    tmp[2] .*= transpose(wnp)

    arr[:, 1:new_m] = tmp[1] + tmp[2]
    arr[:, new_m+1:m] = tmp[1] - tmp[2]
end

function _fft2_pow2_default(arr::Matrix)
    arr = _fft1_pow2(arr, 1)
    arr = _fft1_pow2(arr, 2)
    return arr
end

function _fft2_pow2_cooley(arr::Matrix)
    arr = ComplexF64.(arr)  # make a complex copy
    _fft2_pow2_cooley_rect!(arr)
    return arr
end
