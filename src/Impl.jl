T1 = Union{SubArray{ComplexF64, 1}, Vector{ComplexF64}}
T2 = Union{SubArray{ComplexF64, 2}, Matrix{ComplexF64}}

function intlen(n)
    return UInt64(round(log2(n)))
end

function _rev_k_bits(x::UInt64, k::Int)
    @fastmath begin
        x = ((x >> 32) & 0b0000000000000000000000000000000011111111111111111111111111111111) | ((x & 0b0000000000000000000000000000000011111111111111111111111111111111) << 32)
        x = ((x >> 16) & 0b0000000000000000111111111111111100000000000000001111111111111111) | ((x & 0b0000000000000000111111111111111100000000000000001111111111111111) << 16)
        x = ((x >> 8) & 0b0000000011111111000000001111111100000000111111110000000011111111) | ((x & 0b0000000011111111000000001111111100000000111111110000000011111111) << 8)
        x = ((x >> 4) & 0b0000111100001111000011110000111100001111000011110000111100001111) | ((x & 0b0000111100001111000011110000111100001111000011110000111100001111) << 4)
        x = ((x >> 2) & 0b0011001100110011001100110011001100110011001100110011001100110011) | ((x & 0b0011001100110011001100110011001100110011001100110011001100110011) << 2)
        x = ((x >> 1) & 0b0101010101010101010101010101010101010101010101010101010101010101) | ((x & 0b0101010101010101010101010101010101010101010101010101010101010101) << 1)
        return x >> (64 - k)
    end
end

function _rev_k_bits_and_move(n::UInt64, len::Int, k::Int)
    reversed = _rev_k_bits(n & (2^(k+1) - 1), k)
    return reversed << (len - k) + ((n >> k) & (2^(len-k) - 1))
end

"""
Rearrange 1d-arrays's values in k levels of bit reverse order in-place
"""
function _bit_reverse_permutation_1d!(arr::T1, k::Int)
    n = UInt64(size(arr, 1))
    len = intlen(n)

    for i in 2:n
        j = _rev_k_bits_and_move(i - 1, len, k) + 1
        if i > j
            arr[i], arr[j] = arr[j], arr[i]
        end
    end
end

"""
Rearrange 1d-arrays's values in bit reverse order in-place
"""
function _bit_reverse_permutation_1d!(arr::T1)
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

function _fft1_pow2_process_block2!(arr::T1, block::Int)
    @fastmath @inbounds begin
        arr[block + 1] += arr[block + 2]
        arr[block + 2] = arr[block + 1] - 2*arr[block + 2]
    end
end

function _fft1_pow2_process_block!(arr::T1, half::Int, block::Int, roots::Vector{ComplexF64})
    @fastmath @inbounds @simd for i in 1:half
        arr[block + i + half] *= roots[i]
        arr[block + i] += arr[block + i + half]
        arr[block + i + half] = arr[block + i] - 2*arr[block + i + half]
    end
end

function _fft1_pow2!(arr::T1, roots::Union{Vector{ComplexF64}, Nothing}=nothing)
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
        @fastmath @inbounds @simd for i in 0:2:n-1
            _fft1_pow2_process_block2!(arr, i)
        end
    end
    len = 4
    while len <= n
        half = len >> 1
        @fastmath @inbounds @simd for i in 0:half - 1
            roots[i + 1] = cispi(-2 * i / len)
        end
        @fastmath @inbounds @simd for i in 0:len:n-1
            @fastmath @inbounds _fft1_pow2_process_block!(arr, half, i, roots)
        end
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
function _bit_reverse_permutation_2d!(arr::T2, k::Int, axis::Int)
    m, n = UInt64.(size(arr))
    len = Int(intlen(axis==1 ? n : m))
    arr_copy = copy(arr)

    for i in 1:(axis==1 ? n : m)
        j = _rev_k_bits_and_move(i - 1, len, k) + 1
        for k in 1:(axis==1 ? m : n)
            if axis == 1
                arr[k, j] = arr_copy[k, i]
            else
                arr[j, k] = arr_copy[i, k]
            end
        end
    end
end

"""
Rearrange 2d-arrays's columns and rows in bit reverse order in-place
"""
function _bit_reverse_permutation_2d!(arr::T2, axis::Int=1)
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

function _fft2_pow2_process_block2!(arr::T2, block1::Int, block2::Int)
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
end

function _fft2_pow2_cooley_square!(arr::T2)
    typeof(arr) <: AbstractArray || error("got non-array")
    ndims(arr) == 2              || error("got more than 2-d array")
    ispow2(size(arr, 1))         || ispow2(size(arr, 2)) || error("array's dimensions are not power of 2")

    set_zero_subnormals(true)

    n = size(arr, 1)

    _bit_reverse_permutation_2d!(arr, 1)
    _bit_reverse_permutation_2d!(arr, 2)

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
        @fastmath @inbounds @simd for i in 0:len - 1
            roots[i + 1] = cispi(-2 * i / len)
        end
        @fastmath @inbounds @simd for block2 in 0:len:n-1
            @fastmath @inbounds @simd for j in 1:half
                @fastmath @inbounds @simd for block1 in 0:len:n-1
                    @fastmath @inbounds @simd for i in 1:half
                        @fastmath @inbounds begin
                            u11 = arr[block1 + i, block2 + j]
                            u21 = arr[block1 + i + half, block2 + j] * roots[i]
                            u22 = arr[block1 + i + half, block2 + j + half] * roots[((i + j - 2) & (len - 1)) + 1]
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

"""
More columns
"""
function _fft2_pow2_cooley_rect1!(arr::Matrix{ComplexF64})
    m, n = size(arr)

    levels = Int(intlen(n ÷ m))
    _bit_reverse_permutation_2d!(arr, levels, 1)
    for j in 1:m:n
        @views _fft2_pow2_cooley_square!(arr[:, j:j+m-1])
    end

    len = 2 * m
    while len <= n
        half = len >> 1
        roots = transpose(cispi.(-2 * (0:half-1) / len))
        for block in 0:len:n-1
            for j in 1:half
                @. @views arr[:, block + j + half] *= roots[j]
                @. @views arr[:, block + j] += arr[:, block + j + half]
                @. @views arr[:, block + j + half] = arr[:, block + j] - 2 * arr[:, block + j + half]
            end
        end
        len <<= 1
    end
end

"""
More rows
"""
function _fft2_pow2_cooley_rect2!(arr::Matrix{ComplexF64})
    arr_perm = permutedims(arr, (2, 1))
    _fft2_pow2_cooley_rect1!(arr_perm)
    permutedims!(arr, arr_perm, (2, 1))
end

function _fft2_pow2_default(arr::Matrix)
    arr = _fft1_pow2(arr, 1)
    arr = _fft1_pow2(arr, 2)
    return arr
end

function _fft2_pow2_cooley(arr::Matrix)
    arr = ComplexF64.(arr)  # make a complex copy
    m, n = size(arr)
    if m == n
        _fft2_pow2_cooley_square!(arr)
    elseif n > m
        _fft2_pow2_cooley_rect1!(arr)
    else
        _fft2_pow2_cooley_rect2!(arr)
    end
    return arr
end
