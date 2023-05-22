module Fourier

include("fft_impl.jl")

# Public Interface:

"""
Computes Fast Fourier Transform of 1-d array `arr`
along given axis `axis`
"""
function fft(arr, axis)

end

"""
Computes Fast Fourier Transform of 2-d array `arr`
`algorithm` can be one of:
- "default": use default 1-d FFT along rows than columns
- "cooley": use analog of Cooley-Tuckey algorithm
"""
function fft2(arr; algorithm="default")
    if algorithm == "default"

    elseif algorithm == "cooley"

    end
end

end # module
