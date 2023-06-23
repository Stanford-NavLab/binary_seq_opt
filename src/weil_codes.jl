# WeilCodes
using LinearAlgebra, Statistics, FFTW

function xcor(x::Vector, y::Vector)
    fx = fft(x)
    fy = fft(y)
    return real(ifft(fx .* conj.(fy)))
end

function prime(num)    
    isPrime = true
        for i in 2:convert(Int64, round(num/2, digits = 0))
            if num % i == 0
                isPrime = false
                break
            end
        end
    return isPrime
end

# get legendre sequence (used to create Weil codes)
function get_leg_seq(L::Int)

    # note that L must be a prime number
    if !prime(L)
        error("L must be a prime number")
    end
    
    legendre_set = ones(Int, 1, L-1)
    for i in 1:(L-1)
        val = mod(i*i, L)
        legendre_set[i] = val
    end
    legendre_set = unique(legendre_set)
    
    # First is defined as a (binary) 1 -- i.e., -1 for (+/- 1 regime)
    leg_seq = -1*ones(Int, 1, L)
    leg_seq_flip_i = legendre_set .+ 1
    for i in leg_seq_flip_i
        leg_seq[i] = 1
    end
    
    return leg_seq
    
end

function get_weil_codes(L::Int)
    # returns Weil codes as a set of +/- 1 sequences

    # L must be a prime number for Weil codes
    if !prime(L)
        error("L must be a prime number")
    end
    
    leg_seq = get_leg_seq(L);
    nWeil = Int((L-1)/2);
    weil_codes = zeros(nWeil, L);
    for i = 1:nWeil
        weil_codes[i, :] = leg_seq .* circshift(leg_seq, [0,i]);
    end
    
    return weil_codes

end

# Test function for legendre sequence (check that matches the length-11 and length-31 sequences from Wikipedia)
function test_leg_seq()
    # legendre sequences copy-pasted from Wikipedia (https://en.wikipedia.org/wiki/Legendre_symbol), with 1st number set to be -1
    # (with first one being a -1 instead of "0")
    
    # length 11
    ls_11 = [-1 1 -1 1 1 1 -1 -1 -1 1 -1]
    ls11_test = get_leg_seq(11)
    ls11_matches = (ls_11 == ls11_test)
    
    # length 31
    ls_31 = [-1 1 1 -1 1 1 -1 1 1 1 1 -1 -1 -1 1 -1 1 -1 1 1 1 -1 -1 -1 -1 1 -1 -1 1 -1 -1]
    ls31_test = get_leg_seq(31)
    ls31_matches = (ls_31 == ls31_test)
    
    # final test
    fin_test = (ls11_matches && ls31_matches)    
    println("Test passed? ", fin_test)
    
end


function test_weil_codes(L::Int, tol = 1e-6)
    weil_codebook = get_weil_codes(L);
    n_codes = Int((L-1)/2);

    all_tests_satisfied = true
    # validate provided codebook satisfies size
    if size(weil_codebook, 1) != (n_codes) || size(weil_codebook, 2) != L
        error(
            "weil_codebook size invalid, should be of size (" *
            string(n_codes) *
            ", " *
            string(L) *
            "), but is of size: " *
            string(size(weil_codebook)),
        )
    end

    all_mean_sqr_auto = zeros(1,n_codes);
    max_corr = 0;
    for idx1 = 1:n_codes   # sequence 1 from family
        for idx2 = idx1:n_codes   # sequence 2 from family

            # perform correlation between two sequences
            corr_test = xcor(weil_codebook[idx1, :], weil_codebook[idx2, :]);
            corr_test = round.(corr_test);

            # if indices the same, perform auto-correlation test
            if idx1 == idx2
                # check zero-lag correlation is same as length of Gold
                crct_first = abs(corr_test[1] - L) < tol

                # save mean squared auto-correlation value
                all_mean_sqr_auto[1,idx1] = mean(corr_test[2:L] .* corr_test[2:L]) 

                # check that both zero-lag and side lobe tests satisfied
                curr_satisfied = crct_first 

                # Update whether all tests are satisfied
                all_tests_satisfied = all_tests_satisfied && curr_satisfied

                # save maximum absolute correlation so far
                max_corr = max(max_corr, maximum( abs.(corr_test[2:L]) ) )
            
            else
                # save maximum absolute correlation so far
                max_corr = max(max_corr, maximum( abs.(corr_test) ) )
            end           

        end
    end
    println("Test passed? ", true)
    println("All unique mean-sqr auto-correlations: ", unique(all_mean_sqr_auto))
    println("Maximum absolute correlation: ", max_corr)
    
end

