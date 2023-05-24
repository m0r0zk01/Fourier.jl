using BenchmarkTools
using FFTW
using Fourier2
using Test

sq = rand(4096, 4096)
@btime fft($sq)  # FFTW
@btime fft2_pow2($sq, algorithm="default")  # My default
@btime fft2_pow2($sq, algorithm="cooley")  # My Cooley
