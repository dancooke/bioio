//
//  main.cpp
//
//  Created by Daniel Cooke on 13/10/2014.
//  Copyright (c) 2014 Oxford University. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>

#include "bioio.h"

using std::string;
using std::cout;
using std::endl;
using std::size_t;

void indexed_fasta_example(const string& fasta_path, const string& fasta_index_path)
{
    bioio::FastaIndex the_fasta_index {bioio::read_fasta_index(fasta_index_path)};
    
    std::ifstream fasta {fasta_path};
    
    //auto chr20 = bioio::read_fasta_contig(fasta, the_fasta_index["20"]);
    
    auto chr4_section = bioio::read_fasta_contig(fasta, the_fasta_index["4"], 5e5, 5.1e5);
    
    cout << chr4_section << endl;
}

void fasta_records_example(const string& fasta_path)
{
    std::ifstream fasta {fasta_path};
}

void fastq_records_example(const string& fasta_path)
{
    std::ifstream fastq {fasta_path};
    
    // Read the first 10 records
    auto some_reads = bioio::read_fastq(fastq, 10);
    
    cout << some_reads[0].name << endl;
    cout << some_reads[0].seq << endl;
    cout << some_reads[0].qual << endl;
    
    // Read the rest of the file, but only use the first word in record info line as record name
    auto more_reads = bioio::read_fastq(fastq, [] (const std::string& name) { return name.substr(1, name.find(" ") - 1); });
    
    cout << more_reads[0].name << endl;
    cout << more_reads[0].seq << endl;
    cout << more_reads[0].qual << endl;
}

int main(int argc, const char * argv[])
{
    std::string home_dir {getenv("HOME")};
    std::string genomics_dir {"/Genomics/"};
    std::string reference_dir {"References/"};
    std::string human_reference {"human_g1k_v37.fasta"};
    
    std::string human_reference_fasta {home_dir + genomics_dir + reference_dir + human_reference};
    std::string human_reference_fasta_index {human_reference_fasta + ".fai"};
    
    std::string a_fastq_file {home_dir + genomics_dir + "Nanopore/Data/chr4_2D.fastq"};
    
    indexed_fasta_example(human_reference_fasta, human_reference_fasta_index);
    //fastq_records_example(a_fastq_file);
    
    return 0;
}
