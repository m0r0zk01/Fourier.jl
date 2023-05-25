using BenchmarkTools
using FFTW
using Fourier2
using Test

BenchmarkTools.DEFAULT_PARAMETERS.samples = 10
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 60
BenchmarkTools.DEFAULT_PARAMETERS.gcsample = true

function display(msg, res)
    println(msg)
    io = IOBuffer()
    show(IOContext(io, :color => true), "text/plain", res);
    s = String(take!(io));
    println(s)
    println()
end

function bench1d(n::Int)
    mat = rand(n)
    println("$n vector")
    bench1 = @benchmark fft($mat)
    bench2 = @benchmark fft_pow2($mat)
    display("FFTW:", bench1)
    display("1d Cooley:", bench2)
end

function bench2d(m::Int, n::Int)
    mat = rand(m, n)
    println("$m x $n matrix")
    bench1 = @benchmark fft($mat)
    bench2 = @benchmark fft2_pow2($mat, algorithm="default")
    bench3 = @benchmark fft2_pow2($mat, algorithm="cooley")
    display("FFTW:", bench1)
    display("Default:", bench2)
    display("Cooley:", bench3)
end

# bench1d(16384)
# bench1d(2^25)

bench2d(128, 128)
bench2d(512, 512)
bench2d(1024, 1024)
bench2d(2048, 2048)
bench2d(4096, 4096)
bench2d(8192, 8192)
bench2d(16384, 16384)

bench2d(128, 16384)
bench2d(512, 16384)
bench2d(4096, 16384)
bench2d(16384, 128)
bench2d(16384, 512)
bench2d(16384, 4096)
