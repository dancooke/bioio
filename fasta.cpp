/*  fasta.cpp
 
 Copyright (C) 2015 University of Oxford.
 
 Author: Daniel Cooke <dcooke@well.ox.ac.uk>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.  */

#include <iostream>
#include <fstream>
#include <string>
#include <cstddef>
#include <regex>
#include <stdexcept>
#include <algorithm>
#include <tuple>

#include "bioio.h"

using GenomicRegion = std::tuple<std::string, size_t, size_t>;

GenomicRegion parse_region(const std::string& region, const bioio::FastaIndex& index)
{
    std::string filtered_region {};
    
    std::remove_copy(region.cbegin(), region.cend(), std::back_inserter(filtered_region), ',');
    
    const static std::regex re {"([^:]+)(?::(\\d+)(-)?(\\d*))?"};
    std::smatch match {};
    
    if (std::regex_search(filtered_region, match, re) && match.size() == 5) {
        auto contig_name = match.str(1);
        
        if (index.count(contig_name) == 0) {
            throw std::runtime_error {"contig " + contig_name + " not found"};
        }
        
        auto contig_size = index.at(contig_name).length;
        
        size_t begin {}, end {};
        
        if (match.str(2).empty()) {
            end = contig_size;
        } else {
            begin = static_cast<size_t>(std::stoull(match.str(2)));
            
            if (match.str(3).empty()) {
                end = begin;
            } else if (match.str(4).empty()) {
                end = contig_size;
            } else {
                end = static_cast<size_t>(std::stoull(match.str(4)));
            }
            
            if (begin > contig_size || end > contig_size) {
                throw std::runtime_error {"region " + region + " is larger than contig " + contig_name + ":0-" + std::to_string(contig_size)};
            }
        }
        
        return {std::move(contig_name), begin, end};
    }
    
    throw std::runtime_error {"region" + region + " has invalid format"};
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
        return 0;
    } else if (argc < 3) {
        std::cerr << "Error: not enough command line arguments" << std::endl;
        print_usage();
        return 1;
    }
    
    auto index_path = get_cmd_option(argv, argv + argc, "-i");
    
    std::string fasta_path {argv[argc - 2]};
    
    std::ifstream fasta {fasta_path, std::ios::binary | std::ios::beg};
    
    if (!fasta) {
        std::cerr << "Error: could not open fasta " << fasta_path << std::endl;
        print_usage();
        return 1;
    }
    
    if (index_path.empty()) {
        if (index_path.find(".") != std::string::npos) {
            index_path = fasta_path;
            index_path.replace(index_path.begin() + index_path.find_last_of("."), index_path.end(), ".fai");
        }
    }
    
    std::ifstream index_file {index_path, std::ios::binary | std::ios::beg};
    
    if (!index_file) {
        index_file.open(fasta_path + ".fai");
        if (!index_file) {
            std::cerr << "Error: could not open index file, use samtools faidx <fasta> to make index file" << std::endl;
            return 1;
        }
    }
    
    try {
        auto index = bioio::read_fasta_index(index_file);
        
        auto region        = parse_region(argv[argc - 1], index);
        const auto& contig = std::get<0>(region);
        auto begin         = std::get<1>(region);
        auto length        = std::get<2>(region) - begin;
        
        if (cmd_option_exists(argv, argv + argc, "-s")) {
            std::cout << length << std::endl;
        } else {
            std::cout << bioio::read_fasta_contig(fasta, index.at(contig), begin, length) << std::endl;
        }
        
        return 0;
    } catch (std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}