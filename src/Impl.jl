function _fft1_pow2_inplace(arr)
    typeof(arr) <: AbstractArray || error("got non-array")
    ndims(arr) == 1 || error("got more than 1-d array")
    ispow2(size(arr)[1]) || error("array's length is not power of 2")

    n = size(arr)[1]

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
        for i in 0:len:n-1
            w = 1 + 0im
            for j in 1:half
                u = arr[i + j]
                v = arr[i + j + half] * w
                arr[i + j] = u + v
                arr[i + j + half] = u - v
                w *= wlen
            end
        end
        len <<= 1
    end
end

function _fft1_pow2(arr, axis)
    arr = Complex{Float64}.(arr)  # make a complex copy
    dims = ndims(arr)
    if dims == 1
        _fft1_pow2_inplace(arr)
    else
        idxs = tuple(deleteat!(collect(1:dims), axis)...)
        arr_perm = axis != 1 ? permutedims(arr, (axis, idxs...)) : arr
        slices = eachslice(arr_perm, dims=(2:dims...,))
        _fft1_pow2_inplace.(slices)
        if axis != 1
            permutedims!(arr, arr_perm, (2:axis..., 1, axis+1:dims...))
        end
    end
    return arr
end

function _fft2_pow2_default(arr)
    arr = _fft1_pow2(arr, 1)
    arr = _fft1_pow2(arr, 2)
    return arr
end

function _fft2_pow2_cooley(arr)

end
