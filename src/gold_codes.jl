"""
To get a Gold code family with length L, call:
G = bin_to_pm1(get_gold_code_family(L))'

G is a Lx(L+2) matrix
"""

function get_lfsr_sequence(tap_array::Array{Int}, reg::Array{Int})
    # returns LFSR sequence as binary (0,1) sequence
    #    tap_array -- taps of the register
    #    reg -- initial state of register (can be anything but all 0s)
    #
    # Note: assuming register (reg) will go from right to left
    #    (i.e. index 1 of reg will be first output of register) 
    #    so reg (and tap array) should be initialized accordingly

    # check that inputs are row vectors
    if size(tap_array, 1) != 1 || size(reg, 1) != 1
        error("both tap_array and reg must be row vectors (at least one is not)")
    end

    # check length of taps and initial register
    if length(tap_array) != length(reg)
        error("tap_array and reg not of same length")
    end

    # initialize full maximum length sequence
    reg_len = length(reg)
    seq_len = 2^length(tap_array) - 1
    full_seq = hcat(reg, zeros(Int, 1, seq_len - reg_len))

    # create rest of the LFSR sequence
    for idx = (reg_len+1):seq_len
        # get output of LFSR, put it in the sequence
        nxt_bit = sum(reg .* tap_array) % 2
        full_seq[idx] = nxt_bit

        # update current register
        reg[1:(reg_len-1)] = reg[2:reg_len]
        reg[reg_len] = nxt_bit
    end
    return full_seq
end

function get_gold_code_family_from_taps(tap1::Array{Int}, tap2::Array{Int})
    # returns Gold code family as a codebook of binary (0,1) sequences
    # Note: family is not in any particular order, starting with all relative delays, then adding ML sequences

    # check that everything is a row vector
    if size(tap1, 1) != 1 || size(tap2, 1) != 1
        error("both tap1 and tap2 must be row vectors (at least one is not)")
    end

    # check length of taps
    if length(tap1) != length(tap2)
        error("tap1 and tap2 not of same length")
    end

    # get sequence length and Gold family size
    nbits = 2^length(tap1) - 1
    gold_fam_size = nbits + 2

    # get maximum length sequences
    mls1 = get_lfsr_sequence(tap1, ones(Int, 1, length(tap1)))
    mls2 = get_lfsr_sequence(tap2, ones(Int, 1, length(tap2)))
    mls1_pm1 = bin_to_pm1(mls1)
    mls2_pm1 = bin_to_pm1(mls2)

    # initialize codebook
    gold_codebook = zeros(Int, gold_fam_size, nbits)

    # get relative delay between maximum-length codes
    delays = 1:nbits
    for delay in delays
        # circularly shift MLS 1, perform element-wise multiplication with MLS 2
        circshift_mls1_pm1 = circshift(mls1_pm1, [1, delay])
        cur_gold = mls2_pm1 .* circshift_mls1_pm1

        # save sequence in codebook (as a (0,1) binary sequence)
        gold_codebook[delay, :] = pm1_to_bin(cur_gold)

    end

    # add in original LFSR codes (also as (0,1) binary sequences)
    gold_codebook[nbits+1, :] = mls1
    gold_codebook[nbits+2, :] = mls2

    return gold_codebook
end

function get_mls_pairs(len_gold::Int)
    # Specify preferred pairs (of maximum-length sequence generators) for generating Gold codes
    # For most preferred pairs, following: https://www.gaussianwaves.com/2015/06/gold-code-generator/ 
    # For length-31: got pair from GPS textbook
    # For length-1023: got pair from GPS ICD

    # length-31 sequences
    if len_gold == 31
        tap1 = [1 0 1 0 0]
        tap2 = [1 0 1 1 1]

        # length-63 sequences
    elseif len_gold == 63
        tap1 = [1 0 0 0 0 1]
        tap2 = [1 1 0 0 1 1]

        # length-127 sequences
    elseif len_gold == 127
        tap1 = [1 0 0 0 1 1 1]
        tap2 = [1 0 0 0 1 0 0]

        # length-511 sequences
    elseif len_gold == 511
        tap1 = [1 0 0 1 0 1 1 0 0]
        tap2 = [1 0 0 0 0 1 0 0 0]

        # length-1023 sequences
    elseif len_gold == 1023
        tap1 = [1 0 0 0 0 0 0 1 0 0]
        tap2 = [1 1 1 0 1 0 0 1 1 0]

        # otherwise, if not within specified
    else
        error(
            "cannot create Gold codes for specified length: len_gold =",
            len_gold,
            "\n",
            "(either inconsistent with Gold code definition, or longer than 1023)",
        )
    end
    return tap1, tap2
end

function get_gold_code_family(len_gold::Int)
    # Get MLS pairs
    tap1, tap2 = get_mls_pairs(len_gold)

    # get and return corresponding Gold codes
    gold_codebook = get_gold_code_family_from_taps(tap1, tap2)
    return gold_codebook
end

function test_mls(mls_seq::Array{Int}, len_mls::Int, tol = 1e-6)
    # mls_seq -- (1, len_mls) binary (0,1) sequence
    # len_mls -- must be a number of value 2^N - 1 (do not explicitly check for this)
    if size(mls_seq, 1) != 1 || size(mls_seq, 2) != len_mls
        error(
            "mls_seq must be a row vector of length len_mls (specified to be " *
            string(len_mls) *
            "), but is of size: " *
            string(size(mls_seq)),
        )
    end
    mls_pm1 = bin_to_pm1(mls_seq)
    auto_corr_mls = xcor(vec(mls_pm1), vec(mls_pm1))
    crct_first = abs(auto_corr_mls[1] - len_mls) < tol
    crct_rest = all([abs(auto_corr_mls[i] + 1) < tol for i = 2:len_mls])

    return crct_first && crct_rest
end

function test_all_mls()
    for len_gold in [31, 63, 127, 511, 1023]
        tap1, tap2 = get_mls_pairs(len_gold)
        mls1 = get_lfsr_sequence(tap1, ones(Int, 1, length(tap1)))
        mls2 = get_lfsr_sequence(tap2, ones(Int, 1, length(tap2)))
        chk_mls1 = test_mls(mls1, len_gold)
        chk_mls2 = test_mls(mls2, len_gold)
        println(
            "Tests passed for both MLS sequences of len_gold=" *
            string(length(mls1)) *
            " ? ",
            chk_mls1 && chk_mls2,
        )
    end
end

function unique_gold_corr_vals(reg_len::Int)
    beta_fcn = 1 + 2^(floor((reg_len + 2) / 2))
    return [-1, -1 * beta_fcn, beta_fcn - 2]
end

function test_gold(gold_codebook::Array{Int}, len_reg::Int, tol = 1e-6)
    # gold_codebook -- (K, N) matrix of the binary (0,1) Gold codes 
    #     where K is number of codes and should be N+2 and N is length of codes
    # len_reg -- length of register to create Gold codes (of length 2^(len_reg)-1), 
    #     according to definition of Gold codes len_reg should not be divisible by 4
    #
    # checks that all auto / cross correlations satisfy unique values of Gold codes

    # validate length of register
    if len_reg <= 0 || len_reg % 4 == 0
        error(
            "len_reg invalid, must be positive integer not divisible by 4, but is: " *
            string(len_reg),
        )
    end

    # get length of codes and number of codes
    len_gold = 2^len_reg - 1
    n_codes = len_gold + 2

    # validate provided codebook satisfies size
    if size(gold_codebook, 1) != (n_codes) || size(gold_codebook, 2) != len_gold
        error(
            "gold_codebook size invalid, should be of size (" *
            string(n_codes) *
            ", " *
            string(len_gold) *
            "), but is of size: " *
            string(size(gold_codebook)),
        )
    end

    # get unique set of correlation values for Gold codes with this register length
    corr_vals = unique_gold_corr_vals(len_reg)

    # start with this as true (will be set to false if any 1 criterion not met)
    all_tests_satisfied = true

    # go through each auto / cross correlation
    for idx1 = 1:n_codes   # sequence 1 from family
        for idx2 = idx1:n_codes   # sequence 2 from family

            # perform correlation between two sequences
            corr_test =
                xcor(bin_to_pm1(gold_codebook[idx1, :]), bin_to_pm1(gold_codebook[idx2, :]))
            rdd_corr = round.(corr_test)

            # if indices the same, perform auto-correlation test
            if idx1 == idx2
                # check zero-lag correlation is same as length of Gold
                crct_first = abs(corr_test[1] - len_gold) < tol

                # check that rest of correlations are from the unique Gold correlation set
                crct_rest =
                    all([minimum(abs.(corr_test[i] .- corr_vals)) < tol for i = 2:len_gold])

                # check that both zero-lag and side lobe tests satisfied
                curr_satisfied = crct_first && crct_rest

                # if indices not the same, perform cross-correlation test
            else
                # check that all correlations are from the unique Gold correlation set
                curr_satisfied =
                    all([minimum(abs.(corr_test[i] .- corr_vals)) < tol for i = 1:len_gold])
            end

            # Update whether all tests are satisfied
            all_tests_satisfied = all_tests_satisfied && curr_satisfied

        end
    end

    # return results of all tests
    return all_tests_satisfied
end

function test_all_gold()
    for len_reg in [5, 6, 7, 9, 10]
        code_length = 2^len_reg - 1
        gold_codebook = get_gold_code_family(code_length)
        chk_gold = test_gold(gold_codebook, len_reg)
        println(
            "Tests passed for Gold code family of len_gold=",
            code_length,
            " ? ",
            chk_gold,
        )
    end
end
