module Fourier2

export fft_pow2, fft2_pow2

include("Impl.jl")

# Public Interface:

"""
Computes Fast Fourier Transform of 1-d array `arr` along given axis `axis`

`arr` length must be power of 2
"""
function fft_pow2(arr, axis=0)
    typeof(arr) <: AbstractArray || error("got non-array")

    return _fft1_pow2(arr, axis != 0 ? axis : ndims(arr))
end

"""
Computes Fast Fourier Transform of 2-d array `arr`
`algorithm` can be one of:
- "default": use default 1-d FFT along rows than columns
- "cooley": use analog of Cooley-Tuckey algorithm

`arr` length must be power of 2
"""
function fft2_pow2(arr; algorithm="default")
    typeof(arr) <: AbstractArray || error("got non-array")
    ndims(arr) == 2              || error("arr is not a 2d array")

    if algorithm == "default"
        return _fft2_pow2_default(arr)
    elseif algorithm == "cooley"
        return _fft2_pow2_cooley(arr)
    else
        error("Unknown algorithm name: $algorithm")
    end
end

end # module
