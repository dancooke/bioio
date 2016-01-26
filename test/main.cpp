#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <algorithm>
#include <iterator>

#include "../bioio.hpp"

namespace test_files
{
    static const std::string lambda       {"lambda_ref.fasta"};
    static const std::string lambda_index {"lambda_ref.fasta.fai"};
    static const std::string knucleotide  {"knucleotide.fasta"};
} // namespace paths

bool test_lambda_with_index()
{
    const auto index = bioio::read_fasta_index(test_files::lambda_index);
    
    const std::string contig {"burn-in"};
    
    if (index.size() != 1) return false;
    if (index.count(contig) != 1) return false;
    if (index.at(contig).length != 48502) return false;
    
    const auto sequence = bioio::read_fasta_contig(test_files::lambda, index.at(contig));
    
    if (sequence.size() != 48502) return false;
    
    constexpr size_t chunk_begin {2000};
    constexpr size_t chunk_size  {1000};
    
    const auto chunk = bioio::read_fasta_contig<std::vector<char>>(test_files::lambda, index.at(contig), chunk_begin, chunk_size);
    
    if (chunk.size() != chunk_size) return false;
    if (!std::equal(std::cbegin(chunk), std::cend(chunk), std::next(std::cbegin(sequence), chunk_begin))) return false;
    
    return true;
}

bool test_lambda_without_index()
{
    if (bioio::count_fasta_records(test_files::lambda) != 1) return false;
    
    const auto record = bioio::read_single_contig_fasta(test_files::lambda);
    
    //if (record.contig_name != "burn-in") return false;
    if (record.sequence.size() != 48502) return false;
    
    const auto records = bioio::read_fasta(test_files::lambda);
    
    if (records.size() != 1) return false;
    //if (records.front().contig_name != "burn-in") return false;
    if (records.front().sequence.size() != 48502) return false;
    
    return true;
}

bool test_knucleotide()
{
    return true;
}

int main()
{
    unsigned num_tests_passed {0}, num_tests_failed {0};
    
    if (test_lambda_with_index()) {
        ++num_tests_passed;
    } else {
        std::cout << "Failed 'test_lambda_with_index'" << std::endl;
        ++num_tests_failed;
    }
    
    if (test_lambda_without_index()) {
        ++num_tests_passed;
    } else {
        std::cout << "Failed 'test_lambda_without_index'" << std::endl;
        ++num_tests_failed;
    }
    
    if (test_knucleotide()) {
        ++num_tests_passed;
    } else {
        std::cout << "Failed 'test_knucleotide'" << std::endl;
        ++num_tests_failed;
    }
    
    std::cout << "Passed " << num_tests_passed << " tests" << std::endl;
    std::cout << "Failed " << num_tests_failed << " tests" << std::endl;
    
    return 0;
}