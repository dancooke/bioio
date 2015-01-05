//
//  bioio.h
//
//  Created by Daniel Cooke on 13/10/2014.
//  Copyright (c) 2014 Oxford University. All rights reserved.
//

#ifndef __bioio__bioio__
#define __bioio__bioio__

#include <string>
#include <iostream>
#include <fstream>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <exception>

namespace bioio {

/*=======================================================================================
 TYPEDEFS
 =======================================================================================*/

using ReadId     = std::string;
using DnaString  = std::string;
using QualString = std::string;

template<typename T1=std::string, typename T2=std::string>
struct FastaRecord
{
    T1 id;
    T2 seq;
    
    FastaRecord() = delete;
    template<typename T1_, typename T2_>
    FastaRecord(T1_&& id_, T2_&& seq_) : id(std::forward<T1_>(id_)), seq(std::forward<T2_>(seq_)) {}
};

template<typename T1=std::string, typename T2=std::string, typename T3=std::string>
struct FastqRecord
{
    T1 id;
    T2 seq;
    T3 qual;
    
    FastqRecord() = delete;
    template<typename T1_, typename T2_, typename T3_>
    FastqRecord(T1_&& id_, T2_&& seq_, T3_&& qual_) : id(std::forward<T1_>(id_)), seq(std::forward<T2_>(seq_)), qual(std::forward<T3_>(qual_)) {}
};

template<typename T=std::string>
using ReadIds = std::unordered_set<T>;

template<typename T1=std::string, typename T2=std::string>
using FastaMap = std::unordered_map<T1, T2>;
template<typename T1=std::string, typename T2=std::string, typename T3=std::string>
using FastqMap = std::unordered_map<T1, std::pair<T2, T3>>;

template<typename T1=std::string, typename T2=std::string>
using FastaReads = std::pair<ReadIds<T1>, FastaMap<T1, T2>>;
template<typename T1=std::string, typename T2=std::string, typename T3=std::string>
using FastqReads = std::pair<ReadIds<T1>, FastqMap<T1, T2, T3>>;

/*=======================================================================================
 MISC: functions for getting file information.
 =======================================================================================*/

std::size_t get_num_records(std::ifstream& file, char record_delimiter)
{
    std::size_t num_records = 0;
    std::string first_word;
    while (file >> first_word) {
        if (first_word[0] == record_delimiter) {
            ++num_records;
        }
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    // Uses pre-existing filestream so, need to reset eof flags and rewind.
    file.clear();
    file.seekg(0, std::ios::beg);
    return num_records;
}

/*=======================================================================================
    REFERENCE: Optimised for reading a Fasta file with a *single* record.
 =======================================================================================*/

/*!	@function	Reads a FASTA file with a single record.
	@param  path    FASTA file path.
	@return Just the DNA sequence.
 */
template<typename T=std::string>
T read_reference_seq(std::string path)
{
    std::ifstream fasta(path, std::ios::binary | std::ios::ate);
    std::string id;
    T seq;
    size_t file_len = fasta.tellg();
    fasta.seekg(0, std::ios::beg);
    std::getline(fasta, id, '\n');
    seq.resize(file_len - id.size());
    fasta.read(&seq[0], seq.size());
    auto new_end = std::remove(seq.begin(), seq.end(), '\n');
    seq.resize(new_end - seq.begin() - 1);
    return seq;
}

/*!	@function	Reads a FASTA file with a single record.
	@param  path    FASTA file path.
	@return The header and DNA sequence.
	@note
 */
template<typename T1=std::string, typename T2=std::string>
FastaRecord<T1, T2> read_reference(std::string path)
{
    std::ifstream fasta(path, std::ios::binary | std::ios::ate);
    T1 id;
    T2 seq;
    std::size_t len = fasta.tellg();
    fasta.seekg(0, std::ios::beg);
    std::getline(fasta, id, '\n');
    seq.resize(len - id.size() - 1);
    fasta.read(&seq[0], seq.size());
    std::remove(seq.begin(), seq.end(), '\n');
    return FastaRecord<T1, T2>(id, seq);
}

/*!	@function	Reads a segment of a FASTA file with a single record.
	@param  path    FASTA file path.
    @param  start   Substring start index.
    @param  end     Substring end index.
	@return A substring of a DNA sequence.
 */
template<typename T=std::string>
T read_reference_seq(std::string path, std::size_t start, std::size_t end)
{
    if (start > end) {
        std::cout << "start must be before end." << std::endl;
        return "";
    }
    std::ifstream fasta(path, std::ios::binary | std::ios::beg);
    std::string id;
    T seq;
    std::getline(fasta, id, '\n');
    seq.resize(end - start);
    fasta.seekg(start, std::ios::cur);
    fasta.read(&seq[0], end - start);
    return seq;
}

/*!	@function	Reads a segment of a FASTA file with a single record.
	@param  path    FASTA file path.
	@param  start   Substring start index.
    @param  end     Substring end index.
	@return Header and a substring of a DNA sequence.
 */
template<typename T1=std::string, typename T2=std::string>
FastaRecord<T1, T2> read_reference(std::string path, std::size_t start, std::size_t end)
{
    if (start > end) {
        std::cout << "warning: start > end." << std::endl;
        return std::make_pair("", "");
    }
    std::ifstream fasta(path, std::ios::binary | std::ios::ate);
    T1 id;
    T2 seq;
    std::getline(fasta, id, '\n');
    seq.resize(end - start);
    fasta.seekg(start, std::ios::cur);
    fasta.read(&seq[0], end - start);
    return FastaRecord<T1, T2>(id, seq);
}

/*!	@function	Writes a single fasta record to a file.
	@param  path    FASTA file path.
	@param  data    FASTA record containing header and sequence.
	@return 1 for sucess, 0 for failure.
 */
template<typename T1=std::string, typename T2=std::string>
int write_reference(std::string path, const FastaRecord<T1, T2>& data)
{
    std::ofstream fasta(path, std::ios::out | std::ios::binary);
    if (fasta) {
        fasta << ">" << data.first << "\n";
        fasta << data.second;
        return 0;
    }
    return 1;
}

/*=======================================================================================
 FASTA: For reading multiple record Fasta files.
 =======================================================================================*/

/*!	@function	Reads the next record in a FASTA file.
 @param fasta   An open FASTA file.
 @return    A FastaRecord
 */
template<typename T1=std::string, typename T2=std::string>
FastaRecord<T1, T2> read_fasta_record(std::ifstream& fasta)
{
    T1 id;
    std::getline(fasta, id); // id always a single line
    
    T2 line;
    std::getline(fasta, line);
    
    // The FASTA format is not as simple as FASTQ - the sequence
    // may be broken into multiple lines.
    if (fasta.peek() == '>' || !fasta.good()) {
        return FastaRecord<T1, T2>(id, line);
    } else {
        line.shrink_to_fit();
        std::size_t line_size = line.size();
        
        T2 seq;
        auto it = std::back_inserter(seq);
        std::copy(line.begin(), line.end(), it);
        
        while (fasta.peek() != '>' && fasta.good()) {
            fasta.getline(&line[0], line_size + 1);
            std::copy(line.begin(), line.begin() + fasta.gcount(), it);
        }
        
        return FastaRecord<T1, T2>(id, seq);
    }
}

/*!	@function	Reads Fasta and builds map of read ids after applying given string 
                function to each id.
	@param
	@return
	@note
 */
template<typename T=std::string>
std::vector<T> read_fasta_seqs(std::string path)
{
    std::ifstream fasta(path, std::ios::binary);
    std::vector<T> seqs;
    std::size_t num_records = get_num_records(fasta, '>');
    seqs.reserve(num_records);
    while (num_records) {
        seqs.push_back(read_fasta_record<std::string, T>(fasta).seq);
        --num_records;
    }
    return seqs;
}

/*!	@function	Reads Fasta and builds map of read ids after applying given string 
                function to each id.
	@param
	@return
	@note
 */
template<typename T1=std::string, typename T2=std::string>
std::vector<FastaRecord<T1, T2>> read_fasta_seqs(std::string path)
{
    std::ifstream fasta(path, std::ios::binary);
    std::vector<FastaRecord<T1, T2>> seqs;
    std::size_t num_records = get_num_records(fasta, '>');
    seqs.reserve(num_records);
    while (num_records) {
        seqs.push_back(read_fasta_record<T1, T2>(fasta));
        --num_records;
    }
    return seqs;
}

/*!	@function	Reads Fasta and builds map of read ids after applying given string 
                function to each id.
	@param
	@return
	@note
 */
template <typename T1=std::string, typename T2=std::string, typename F>
FastaReads<T1, T2> read_fasta_map(std::string path, F f)
{
   std::ifstream fasta(path, std::ios::binary);
   FastaMap<T1, T2> records;
   ReadIds<T1> ids;
   std::size_t num_records = get_num_records(fasta, '>');
   records.reserve(num_records);
   ids.reserve(num_records);
   while (num_records) {
       auto record = read_fasta_record<T1, T2>(fasta);
       auto f_id = f(std::move(record.id));
       records.emplace(f_id, std::move(record.seq));
       ids.insert(std::move(f_id));
       --num_records;
   }
   return std::make_pair(ids, records);
}

/*!	@function	Reads
	@param
	@return
	@note
 */
template<typename T1=std::string, typename T2=std::string>
inline FastaReads<T1, T2> read_fasta_map(std::string path)
{
    return read_fasta_map(path, [] (ReadId id) { return id; });
}

/*!	@function	Searches a Fasta file for the given record ids, after applying a string function
                to each record id found in the Fasta.
	@param
	@return
	@note
 */
template <typename T1=std::string, typename T2=std::string, typename F>
FastaReads<T1, T2> read_fasta_map(std::string path, const ReadIds<T1>& ids, F f)
{
   std::ifstream fasta(path, std::ios::binary);
   FastaMap<T1, T2> records;
   ReadIds<T1> f_ids;
   std::size_t num_records = ids.size();
   records.reserve(num_records);
   f_ids.reserve(num_records);
   while (num_records) {
       auto record = read_fasta_record<T1, T2>(fasta);
       auto f_id = f(std::move(record.id));
       if (ids.find(f_id) != ids.end()) {
           records.emplace(f_id, std::move(record.seq));
           f_ids.insert(std::move(f_id));
           --num_records;
       }
   }
   return std::make_pair(f_ids, records);
}

/*!	@function	Searches a Fasta file for the given record ids.
	@param
	@return
	@note
 */
template<typename T1=std::string, typename T2=std::string>
inline FastaReads<T1, T2> read_fasta_map(std::string path, const ReadIds<T1>& ids)
{
    return read_fasta_map(path, ids, [] (ReadId id) { return id; });
}

/*!	@function	Reads
	@param
	@return
	@note
 */
template<typename T1=std::string, typename T2=std::string>
int write_fasta(std::string path, const FastaReads<T1, T2>& data)
{
    std::ofstream fasta(path, std::ios::out | std::ios::binary);
    if (fasta) {
        for (auto id : data.first) {
            fasta << ">" << id << "\n";
            fasta << data.second.at(id);
        }
        return 0;
    }
    return 1;
}

/*=======================================================================================
 FASTQ: For reading multiple line Fastq files.
 =======================================================================================*/

/*!	@function	Reads
 @param
 @return
 @note
 */
template<typename T1=std::string, typename T2=std::string, typename T3=std::string>
FastqRecord<T1, T2, T3> read_fastq_record(std::ifstream& fastq)
{
    T1 id;
    T2 seq;
    T3 quals;
    // Unlike FASTA, FASTQ always use one line per field.
    std::getline(fastq, id);
    std::getline(fastq, seq);
    fastq.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::getline(fastq, quals);
    return FastqRecord<T1, T2, T3>(std::move(id), std::move(seq), std::move(quals));
}

/*!	@function	Reads
	@param
	@return
	@note
 */
template<typename T=std::string>
std::vector<T> read_fastq_seqs(std::string path)
{
    std::ifstream fastq(path, std::ios::binary);
    std::vector<T> seqs;
    std::size_t num_records = get_num_records(fastq, '@');
    seqs.reserve(num_records);
    while (num_records) {
        seqs.push_back(read_fastq_record<std::string, T, std::string>(fastq).seq);
        --num_records;
    }
    return seqs;
}

/*!	@function	Reads
	@param
	@return
	@note
 */
template<typename T1=std::string, typename T2=std::string, typename T3=std::string>
std::vector<FastqRecord<T1, T2, T3>> read_fastq(std::string path)
{
   std::ifstream fastq(path, std::ios::binary);
   std::vector<FastqRecord<T1, T2, T3>> data;
   std::size_t num_records = get_num_records(fastq, '@');
   data.reserve(num_records);
   while (num_records) {
       data.push_back(read_fastq_record<T1, T2, T3>(fastq));
       --num_records;
   }
   return data;
}

/*!	@function	Reads
	@param
	@return
	@note
 */
template<typename T1=std::string, typename T2=std::string, typename T3=std::string>
FastqReads<T1, T2, T3> read_fastq_map(std::string path)
{
   std::ifstream fastq(path, std::ios::binary);
   ReadIds<T1> ids;
   FastqMap<T1, T2, T3> data;
   std::size_t num_records = get_num_records(fastq, '@');
   data.reserve(num_records);
   while (num_records) {
       auto record = read_fastq_record<T1, T2, T3>(fastq);
       data.emplace(record.id, std::make_pair(std::move(record.seq), std::move(record.qual)));
       ids.insert(std::move(record.id));
       --num_records;
   }
   return std::make_pair(ids, data);
}

}

#endif /* defined(__bioio__bioio__) */
