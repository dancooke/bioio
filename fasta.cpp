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
        size_t begin {}, end {};
        
        auto contig_size = index.at(contig_name).length;
        
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
                throw std::runtime_error {"region " + region + " is larger than contig"};
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
    if (argc < 3) {
        print_usage();
        return 1;
    }
    
    std::string index_path {};
    
    index_path = get_cmd_option(argv, argv + argc, "-i");
    
    std::string fasta_path {argv[argc - 2]};
    
    if (index_path.empty()) {
        index_path = fasta_path;
        index_path.replace(index_path.begin() + index_path.find_last_of("."), index_path.end(), ".fai");
    }
    
    std::ifstream index_file {index_path};
    
    if (!index_file) {
        index_file.open(fasta_path + ".fai");
        if (!index_file) {
            std::cerr << "Could not open index file" << std::endl;
            return 1;
        }
    }
    
    auto index = bioio::read_fasta_index(index_file);
    
    auto region = parse_region(argv[argc - 1], index);
    
    const auto& contig = std::get<0>(region);
    auto begin         = std::get<1>(region);
    auto length        = std::get<2>(region) - begin;
    
    auto sequence = bioio::read_fasta_contig(fasta_path, index.at(contig), begin, length);
    
    if (cmd_option_exists(argv, argv + argc, "-s")) {
        std::cout << sequence.size() << std::endl;
    } else {
        std::cout << sequence << std::endl;
    }
    
    return 0;
}
