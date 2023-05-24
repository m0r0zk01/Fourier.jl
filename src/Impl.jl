function _fft1_pow2!(arr)
    typeof(arr) <: AbstractArray || error("got non-array")
    ndims(arr) == 1 || error("got more than 1-d array")
    ispow2(size(arr, 1)) || error("array's length is not power of 2")

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

    len = 2
    while len <= n
        half = len >> 1
        ang = -2 * pi / len
        wlen = cos(ang) + 1im * sin(ang)
        for block in 0:len:n-1
            w = 1 + 0im
            for i in 1:half
                u = arr[block + i]
                v = arr[block + i + half] * w
                arr[block + i] = u + v
                arr[block + i + half] = u - v
                w *= wlen
            end
        end
        len <<= 1
    end

    return nothing
end

function _fft1_pow2(arr, axis)
    arr = Complex{Float64}.(arr)  # make a complex copy
    dims = ndims(arr)
    if dims == 1
        _fft1_pow2!(arr)
    else
        idxs = tuple(deleteat!(collect(1:dims), axis)...)
        arr_perm = axis != 1 ? permutedims(arr, (axis, idxs...)) : arr
        slices = eachslice(arr_perm, dims=(2:dims...,))
        _fft1_pow2!.(slices)
        if axis != 1
            permutedims!(arr, arr_perm, (2:axis..., 1, axis+1:dims...))
        end
    end
    return arr
end

function _bit_reverse_permutation_2d!(arr)
    m, n = size(arr)
    j = 1
    for i in 2:n
        bit = n >> 1
        while j > bit
            j -= bit
            bit >>= 1
        end
        j += bit
        if i < j
            for k in 1:m
                arr[k, i], arr[k, j] = arr[k, j], arr[k, i]
            end
        end
    end
end

function _fft2_pow2_cooley_square!(arr)
    typeof(arr) <: AbstractArray || error("got non-array")
    ndims(arr) == 2 || error("got more than 2-d array")
    ispow2(size(arr, 1)) || ispow2(size(arr, 2)) || error("array's dimensions are not power of 2")

    n = size(arr, 1)

    _bit_reverse_permutation_2d!(arr)
    arr_perm = permutedims(arr, (2, 1))
    _bit_reverse_permutation_2d!(arr_perm)
    permutedims!(arr, arr_perm, (2, 1))

    len = 2
    while len <= n
        half = len >> 1
        ang = -2 * pi / len
        wlen = cos(ang) + 1im * sin(ang)
        for block2 in 0:len:n-1
            for block1 in 0:len:n-1
                w1 = 1 + 0im
                w12 = 1 + 0im
                for j in 1:half
                    w2 = 1 + 0im
                    for i in 1:half
                        u11 = arr[block1 + i, block2 + j]
                        u12 = arr[block1 + i, block2 + j + half] * w1
                        u21 = arr[block1 + i + half, block2 + j] * w2
                        u22 = arr[block1 + i + half, block2 + j + half] * w1 * w2

                        arr[block1 + i, block2 + j] = u11 + u12 + u21 + u22
                        arr[block1 + i, block2 + j + half] = u11 - u12 + u21 - u22
                        arr[block1 + i + half, block2 + j] = u11 + u12 - u21 - u22
                        arr[block1 + i + half, block2 + j + half] = u11 - u12 - u21 + u22

                        w2 *= wlen
                        w12 *= wlen
                    end
                    w1 *= wlen
                    w12 *= wlen
                end
            end
        end
        len <<= 1
    end
end

function W(n, p)
    return exp(-2 * 1im * π * p / n)
end

function _fft2_pow2_cooley_square_req!(arr)
    n = size(arr, 1)
    if n <= 2
        if n == 1
            return
        end
        tmp = [
            arr[1, 1] + arr[1, 2] + arr[2, 1] + arr[2, 2] arr[1, 1] - arr[1, 2] + arr[2, 1] - arr[2, 2];
            arr[1, 1] + arr[1, 2] - arr[2, 1] - arr[2, 2] arr[1, 1] - arr[1, 2] - arr[2, 1] + arr[2, 2]
        ]
        copy!(arr, tmp)
        return
    end

    new_n = n >> 1
    tmp = fill(Matrix{ComplexF64}(undef, new_n, new_n), 4)
    tmp[1] = arr[1:2:n, 1:2:n]
    tmp[2] = arr[2:2:n, 1:2:n]
    tmp[3] = arr[1:2:n, 2:2:n]
    tmp[4] = arr[2:2:n, 2:2:n]

    _fft2_pow2_cooley_square_req!.(tmp[1:4])

    println.(tmp)

    wnp = W.(n, 0:new_n-1)
    tmp[2] .*= wnp
    tmp[3] .*= transpose(wnp)
    tmp[4] .*= transpose(wnp)
    tmp[4] .*= wnp

    arr[1:new_n, 1:new_n] = tmp[1] + tmp[2] + tmp[3] + tmp[4]

    arr[new_n+1:n, 1:new_n] = tmp[1] - tmp[2] + tmp[3] - tmp[4]

    arr[1:new_n, new_n+1:n] = tmp[1] + tmp[2] - tmp[3] - tmp[4]

    arr[new_n+1:n, new_n+1:n] = tmp[1] - tmp[2] - tmp[3] + tmp[4]

    nothing
end

function _fft2_pow2_default(arr)
    arr = _fft1_pow2(arr, 1)
    arr = _fft1_pow2(arr, 2)
    return arr
end

function _fft2_pow2_cooley(arr)
    arr = Complex{Float64}.(arr)  # make a complex copy
    _fft2_pow2_cooley_square!(arr)
    return arr
end
