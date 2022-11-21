
""" cross correlation in frequency domain """
function xcor(x::Vector, y::Vector)
    fx = fft(x)
    fy = fft(y)
    return real(ifft(fx .* conj.(fy)))
end

function xcor(x::Vector{Int}, y::Vector{Int})::Vector{Int}
    fx = fft(x)
    fy = fft(y)
    return convert(Vector{Int}, round.(real(ifft(fx .* conj.(fy)))))
end


""" autocorrelation """
function acor(x::Vector)
    fx = fft(x)
    return real(ifft(fx .* conj.(fx)))
end

function acor(x::Vector{Int})::Vector{Int}
    fx = fft(x)
    return convert(Vector{Int}, round.(real(ifft(fx .* conj.(fx)))))
end

""" generate a random Bern(1/2) binary matrix/vector """
function randb(size...)
    return rand((-1, 1), size)
end


""" convert 0/1 arrays to -1/+1 arrays """
function bin_to_pm1(x)
    return -2 * x .+ 1
end


""" convert -1/+1 arrays to 0/1 arrays """
function pm1_to_bin(x)
    return ((x .+ 1) .÷ -2) .+ 1
end
