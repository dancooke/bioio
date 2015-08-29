# Bioio

Bioio is a small, fast, C++11 library for FASTA/Q I/O.

Template methods allow any contiguous containers (e.g. std::vector and std::string).

fasta.cpp is a simple command line application that uses bioio to fetch subsequences from a fasta file. It can be compiled with `make`, and is used like

    ./fasta [options] <fasta_path> <region>

where `options` can include `-i` to set the index path, `-s` to output the length of the given sequence.

`region` is of the format `contig:[begin]-[end]` where `begin` and `end` are optional (e.g. `X` is all of contig X, and `X:1000-` is everything from 1000 to the end of X).

e.g.

    ./fasta human_g1k_v37.fasta X:1,000,000-1,000,100
    ./fasta human_g1k_v37.fasta 5:10,040,010-
