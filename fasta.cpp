/*  fasta.cpp

 Copyright (C) 2017 University of Oxford.

 Author: Daniel Cooke <dcooke@well.ox.ac.uk>

 Use of this source code is governed by the MIT license that can be found in the LICENSE file  */

#include <iostream>
#include <fstream>
#include <string>
#include <cstddef>
#include <regex>
#include <stdexcept>
#include <algorithm>
#include <tuple>

#include "bioio.hpp"

using GenomicRegion = std::tuple<std::string, size_t, size_t>;

GenomicRegion parse_region(std::string region, const bioio::FastaIndex& index)
{
    region.erase(std::remove(region.begin(), region.end(), ','), region.end());

    const static std::regex re {"([^:]+)(?::(\\d+)(-)?(\\d*))?"};

    std::smatch match;

    if (std::regex_match(region, match, re) && match.size() == 5) {
        auto contig_name = match.str(1);

        if (index.count(contig_name) == 0) {
            throw std::runtime_error {"contig " + contig_name + " not found"};
        }

        const auto contig_size = index.at(contig_name).length;

        size_t begin {0}, end {0};

        if (match.str(2).empty()) {
            end = contig_size;
        } else {
            begin = static_cast<size_t>(std::stoull(match.str(2)));

            if (match.str(3).empty()) {
                end = begin + 1;
            } else if (match.str(4).empty()) {
                end = contig_size;
            } else {
                end = static_cast<size_t>(std::stoull(match.str(4)));
            }

            if (begin > contig_size) {
                throw std::runtime_error {"region " + region + " is larger than contig " + contig_name + ":0-" + std::to_string(contig_size)};
            }

            if (begin > end) {
                throw std::runtime_error {"begin position is past end position in region " + region};
            }

            if (end > contig_size) {
                end = contig_size;
            }
        }

        return GenomicRegion {std::move(contig_name), begin, end};
    }

    throw std::runtime_error {"could not parse region " + region};
}

bool cmd_option_exists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

std::string get_cmd_option(char** begin, char** end, const std::string& option)
{
    auto it = std::find(begin, end, option);
    if (it != end && ++it != end) return *it;
    return {};
}

void print_usage()
{
    std::cout << "Usage: fasta [options] <path> <region>" << std::endl;
}

int main(int argc, char **argv)
{
    if (cmd_option_exists(argv, argv + argc, "-h")) {
        print_usage();
        std::cout << "options:" << '\n';
        std::cout << "\t-h\tprint help" << '\n';
        std::cout << "\t-s\toutput sequence size of region" << '\n';
        //std::cout << "\t" << '\n';
        return 0;
    } else if (argc < 3) {
        std::cerr << "Error: not enough command line arguments" << std::endl;
        print_usage();
        return 1;
    }

    auto index_path = get_cmd_option(argv, argv + argc, "-i");

    std::string fasta_path {argv[argc - 2]};

    std::ifstream fasta {fasta_path, std::ios::binary};

    if (!fasta) {
        std::cerr << "Error: could not open fasta " << fasta_path << std::endl;
        print_usage();
        return 1;
    }

    if (index_path.empty() && index_path.find(".") != std::string::npos) {
        index_path = fasta_path;
        index_path.replace(index_path.begin() + index_path.find_last_of("."), index_path.end(), ".fai");
    }

    std::ifstream index_file {index_path, std::ios::binary};

    if (!index_file) {
        index_file.open(fasta_path + ".fai");
        if (!index_file) {
            std::cerr << "Error: could not open index file, use samtools faidx <fasta> to make index file" << std::endl;
            return 1;
        }
    }

    try {
        const auto index = bioio::read_fasta_index(index_file);

        const auto region  = parse_region(argv[argc - 1], index);
        const auto& contig = std::get<0>(region);
        const auto begin   = std::get<1>(region);
        const auto length  = std::get<2>(region) - begin;

        if (cmd_option_exists(argv, argv + argc, "-s")) {
            std::cout << length << std::endl;
        } else {
            const auto sequence = bioio::read_fasta_contig(fasta, index.at(contig), begin, length);

            std::cout << sequence << std::endl;
        }

        return 0;
    } catch (std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
