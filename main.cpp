//
//  main.cpp
//
//  Created by Daniel Cooke on 13/10/2014.
//  Copyright (c) 2014 Oxford University. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <cstdio>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "bioio.h"

using std::string;
using std::cout;
using std::endl;
using std::size_t;

/*=======================================================================================
READ tests: a number of ways to read a FASTA file in C++
 =======================================================================================*/

void test_it(string ref_path)
{
    auto start = std::chrono::system_clock::now();
    
    std::ifstream fasta(ref_path, std::ios::ate);
    
    string seq, header;
    size_t file_len = fasta.tellg();
    fasta.seekg(0, std::ios::beg);
    
    std::getline(fasta, header, '\n');
    size_t seq_len = file_len - header.length();
    seq.reserve(seq_len);
    
    string line;
    while (fasta) {
        std::getline(fasta, line, '\n');
        seq += line;
        fasta.peek();
    }
    
    seq.shrink_to_fit();
    
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    //cout << seq.length() << endl;
    
    cout << "it: " << elapsed.count() << endl;
}

void test_get(string ref_path)
{
    auto start = std::chrono::system_clock::now();
    
    std::ifstream fasta(ref_path, std::ios::ate);
    
    string seq, header;
    
    size_t file_len = fasta.tellg();
    fasta.seekg(0, std::ios::beg);
    
    std::getline(fasta, header, '\n');
    size_t seq_len = file_len - header.length();
    seq.reserve(seq_len);
    
    char c = fasta.get();
    while (fasta) {
        if (c != '\n') {
            seq.push_back(c);
        }
        c = fasta.get();
    }
    
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    //cout << seq.length() << endl;
    
    cout << "cpp_get: " << elapsed.count() << endl;
}

void test_c_get(string ref_path)
{
    auto start = std::chrono::system_clock::now();
    
    std::FILE* fasta = std::fopen(ref_path.c_str(), "r");
    
    string seq, header;
    std::fseek(fasta, 0, SEEK_END);
    size_t file_len = std::ftell(fasta);
    std::rewind(fasta);
    
    while (char c = std::getc(fasta) != '\n') {
        header.push_back(c);
    }
    
    size_t seq_len = file_len - header.length();
    seq.reserve(seq_len);
    
    while (char c = std::getc(fasta) != EOF) {
        if (c != '\n') {
            seq.push_back(c);
        }
    }
    
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    //cout << seq.length() << endl;
    
    cout << "c_get: " << elapsed.count() << endl;
}

void test_remove(string ref_path)
{
    auto start = std::chrono::system_clock::now();
    
    std::ifstream fasta(ref_path, std::ios::binary | std::ios::ate);
    
    std::string seq, header;
    size_t file_len = fasta.tellg();
    fasta.seekg(0, std::ios::beg);
    std::getline(fasta, header, '\n');
    seq.resize(file_len - header.length());
    
    fasta.read(&seq[0], seq.length());
    auto new_end = std::remove(seq.begin(), seq.end(), '\n');
    
    seq.resize(new_end - seq.begin() - 1);
    
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    //cout << seq.length() << endl;
    
    cout << "remove: " << elapsed.count() << endl;
}

void test_maybe_remove(string ref_path)
{
    auto start = std::chrono::system_clock::now();
    
    std::ifstream fasta(ref_path, std::ios::binary | std::ios::ate);
    
    std::string seq, header;
    size_t file_len = fasta.tellg();
    fasta.seekg(0, std::ios::beg);
    std::getline(fasta, header, '\n');
    seq.resize(file_len - header.length() - 1);
    
    fasta.read(&seq[0], seq.length());
    
    auto test_reg_end = seq.cbegin() + 200;
    if (std::find(seq.cbegin(), test_reg_end, '\n') != test_reg_end) {
        auto new_end = std::remove(seq.begin(), seq.end(), '\n');
        seq.resize(new_end - seq.begin());
    }
    
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    //cout << seq.length() << endl;
    
    cout << "maybe remove: " << elapsed.count() << endl;
}

void test_buff(string ref_path)
{
    auto start = std::chrono::system_clock::now();
    
    std::ifstream fasta(ref_path, std::ios::binary | std::ios::ate);
    
    std::string seq, header;
    size_t file_len = fasta.tellg();
    fasta.seekg(0, std::ios::beg);
    std::getline(fasta, header, '\n');
    size_t seq_len = file_len - header.length() - 1;
    seq.reserve(seq_len);
    
    string line;
    std::getline(fasta, line, '\n');
    line.shrink_to_fit();
    size_t line_size = line.length();
    
    auto it = std::back_inserter(seq);
    std::copy(line.begin(), line.end(), it);
    
    if (line_size != seq_len) {
        size_t num_chars_read {};
        
        while (fasta) {
            fasta.read(&line[0], line_size);
            num_chars_read = fasta.gcount();
            std::copy(line.begin(), line.end(), it);
            fasta.ignore();
        }
        
        size_t num_bad_bits = line_size - num_chars_read + 1;
        seq.erase(seq.end() - num_bad_bits - 1, seq.end());
    }
    
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    //cout << seq.length() << endl;
    
    cout << "buff: " << elapsed.count() << endl;
}

void test_bioio(string ref_path)
{
    std::ifstream fasta(ref_path);
    auto index = bioio::read_fasta_index(ref_path + ".fai")["5"];
    
    auto start = std::chrono::system_clock::now();
    
    auto seq = bioio::read_fasta_contig(fasta, index, 100000, 100);
    
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    
    cout << "bioio: " << elapsed.count() << endl;
}

void test_seqan(string ref_path)
{
    auto start = std::chrono::system_clock::now();
    
    seqan::CharString id;
    seqan::Dna5String seq;
    
    seqan::SequenceStream seqStream(ref_path.c_str());
    readRecord(id, seq, seqStream);
    cout << id << '\t' << length(seq) << endl;
    
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "seqan: " << elapsed.count() << endl;
}

void performance_tests()
{
    std::string homedir {getenv("HOME")};
    string ref_path {homedir + "/Genomics/References/human_g1k_v37.fasta"};
    //string ref_path {"/Users/dcooke/Genomics/References/R00000042.fasta"};
    //string ref_path {"/Users/dcooke/Genomics/References/lambda_ref.fasta"};
    
    // test_it(ref_path);
    //
    // test_get(ref_path);
    //
    // test_c_get(ref_path);
    //
    // test_remove(ref_path);
    //
    // test_maybe_remove(ref_path);
    //
    // test_buff(ref_path);
    
    test_bioio(ref_path);
    
    //test_seqan(ref_path);
}

/*=======================================================================================
 Examples: some examples of using bioio
 =======================================================================================*/

void test_read_fasta(string ref_path,  string fasta_path)
{   
//    auto ref  = bioio::read_reference_seq<>(ref_path); // default std::string
//    auto ref1 = bioio::read_reference_seq<std::vector<char>>(ref_path);
//    auto ref2 = bioio::read_reference_seq<>(ref_path, 10000, 20000);
//    
//    auto reads  = bioio::read_fasta<>(fasta_path);
//    
//    auto read_map  = bioio::read_fasta_map<>(fasta_path);
//    // Can pass lambda to strip discription off read id.
//    auto read_map2 = bioio::read_fasta_map<>(fasta_path, [] (std::string id) { return id.substr(0, id.find(' ')); });
    
    // for (auto& id : read_map2.first) {
    //     cout << id << '\n' << read_map2.second.at(id) << endl;
    // }
}

void test_read_fastq(string fastq_path)
{
    // auto reads = bioio::read_fastq_map<>(fastq_path);
    // for (auto& id : reads.first) {
    //     cout << id << '\n' << reads.second.at(id).first << endl;
    // }
}

int main(int argc, const char * argv[])
{   
    performance_tests();
    
//    string ref_path {"/Users/dcooke/Genomics/References/R00000042.fasta"};
//    string fasta_path {"/Users/dcooke/Genomics/Nanopore/Data/tranche1/processed/m1/fastq/R7_tranche1_m1_template.fa"};
//    string fastq_path {"/Users/dcooke/Genomics/Nanopore/Data/lambdaR7/m1/raw_reads/template.fq"};
//    
//    test_read_fasta(ref_path, fasta_path);
//    test_read_fastq(fastq_path);
    
    return 0;
}
