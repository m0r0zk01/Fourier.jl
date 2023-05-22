function _fft1d(arr)
    println("kek")
    typeof(arr) <: AbstractArray || error("_fft1d got non-array")
    ndims(arr) == 1 || error("_fft1d got more than 1-d array")

    n = size(arr)[1]
    arr = convert(Array{Complex}, copy(arr))

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
        half = len ÷ 2
        ang = -2 * pi / len
        wlen = cos(ang) + 1im * sin(ang)
        i = 0
        while i < n
            w = 1 + 0im
            for j in 1:half
                u = arr[i + j]
                v = arr[i + j + half] * w
                arr[i + j] = u + v
                arr[i + j + half] = u - v
                w *= wlen
            end
            i += len
        end
        len <<= 1
    end

    return arr[1:n]
end

function _fft1(arr, axis=1)
    n = size(arr, axis)
    n = nextpow(2, n)

end

function _fft2_default(arr)

end

function _fft2_cooley(arr)

end
