# binary_seq_opt
Binary sequence set optimization for CDMA applications.

Example scripts may be found in the ```icassp23``` and ```sherlock_demo``` directories.

For example, running
```
julia icassp23.jl 0 1023 4 20 ISL 1000
```
finds a binary sequence set with 4 sequences of length 1023. The arguments are, in order, 1) the random seed used to generate the initial sequence set, 2) the sequence lengths, 3) the number of sequences, 4) the block coordinate descent variable subset size, 5) the objective function, and 5) the maximum number of iterations. Additional arguments may be provided, as specified in the scripts.
