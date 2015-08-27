# Bioio

Bioio is a small, fast, C++11 library for FASTA and FASTQ IO.

Template methods allow IO into any continuous containers (e.g. std::vector and std::string).

example.cpp shows some basic use cases of the API.

fasta.cpp is a simple command line application that uses bioio to fetch subsequences from a fasta file. It can be compiled with `make`, and is used like

    ./fasta [options] <fasta_path> <region>

where `options` can include `-i` to set the index path, `-s` to output the length of the given sequence.

`region` is of the format `contig:[begin]-[end]` where `begin` and `end` are optional (e.g. `X` is all of contig X, and `X:1000-` is everything from 1000 to the end of X).
