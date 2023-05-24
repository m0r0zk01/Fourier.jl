"""
(-p)-th n-root of unity
"""
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
