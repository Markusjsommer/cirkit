//
// Created by marku on 10/27/2021.
//

#include "build_index.h"
#include <sstream>
#include <fstream>

//#include <phmap_dump.h>
//#include <parallel_hashmap/phmap_dump.h>

#include <cereal/archives/binary.hpp>
#include "cereal/types/unordered_map.hpp"
#include "cereal/types/memory.hpp"
#include "cereal/types/bitset.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/tuple.hpp"

#include <chrono>
#include <functional>
#include <cstdio>
#include <cinttypes>
#include <iostream>
#include <bitset>
#include <algorithm>
#include "BitsetManager.h"

//template <typename T> using milliseconds = std::chrono::duration<T, std::milli>;

void build_index::read_ref(std::string &gtf, std::string &ref) {
    // read fasta
    std::cout << "Reading fasta at " << ref << std::endl;
    _io.read_fasta(ref, _seq_vec, _contigname_vec);

    // capitalize all nucleotides
    for (auto &seq : _seq_vec){
        for(auto &c: seq){
            c = toupper(c);
        }
    }
    if (_seq_vec.empty()) {
        std::cout << "Fasta read failed..." << std::endl;
        return;
    }

    // read gtf
    std::cout << "Reading gtf at " << gtf << std::endl;
    _gtfio.read_gtf(gtf, _gtf_lines);
}

//void build_index::extract_features(const std::string& feature, bool noalt) {
//    // split gtf lines on tabs
//    // gtf format
//    // 0: contig, 1: source, 2: feature, 3: start, 4: end, 5: score, 6: strand, 7: frame, 8: geneID
//
//    for (auto & l : _gtf_lines){
//        std::istringstream ss(l);
//        std::string substring;
//        std::tuple<std::string,
//                   std::string,
//                   std::string,
//                   std::uint32_t,
//                   std::uint32_t,
//                   std::string,
//                   std::string,
//                   std::string,
//                   std::string> row_tuple =  std::make_tuple("", "", "", 0, 0, "", "", "", "");
//        int i = 0;
//        while (getline(ss, substring, '\t')){
//            // TODO replace with case switch or just avoid this entirely
//            if (i==0){
//                std::get<0>(row_tuple) = substring;
//            } else if (i==1) {
//                std::get<1>(row_tuple) = substring;
//            } else if (i==2) {
//                std::get<2>(row_tuple) = substring;
//            } else if (i==3) {
//                uint32_t x = std::stoi(substring);
//                std::get<3>(row_tuple) = x;
//            } else if (i==4) {
//                uint32_t x = std::stoi(substring);
//                std::get<4>(row_tuple) = x;
//            } else if (i==5) {
//                std::get<5>(row_tuple) = substring;
//            } else if (i==6) {
//                std::get<6>(row_tuple) = substring;
//            } else if (i==7) {
//                std::get<7>(row_tuple) = substring;
//            } else if (i==8) {
//                std::get<8>(row_tuple) = substring;
//            }
//            i++;
//        }
//
//        // skip alt scaffolds from gtf if requested
//        if (noalt){
//            std::string contig;
//            bool alt = false;
//            contig = std::get<0>(row_tuple);
//            if (contig.find("alt") != std::string::npos){
//                alt = true;
//            } else if (contig.find('_') != std::string::npos) {
//                alt = true;
//            }
//            if (alt){
//                continue;
//            }
//        }
//
//        // only keep rows with feature type matching given feature
//        if (std::get<2>(row_tuple) == feature){
//            _gtf_data.emplace_back(row_tuple);
//        }
//    }
//
//    // get sequence on both strands for each feature from the correct contig
//    std::string contigname;
//    int contigname_idx;
//    uint32_t start;
//    uint32_t end;
//    std::string seq_forward;
//    std::string seq_forwardcomplement;
//    std::string geneID;
//    for (auto & x: _gtf_data){
//        contigname = std::get<0>(x);
//        start = std::get<3>(x);
//        end = std::get<4>(x);
//
//        auto it = std::find(_contigname_vec.begin(), _contigname_vec.end(), contigname);
//        if (it == _contigname_vec.end()){
//            std::cout << "contig " << contigname << " from annotation not found in reference fasta, consider removing alt scaffolds with --no-alt" << std::endl;
//            continue;
//        }
//
//        contigname_idx = it - _contigname_vec.begin();
//        seq_forward = _seq_vec[contigname_idx].substr(start, end - start + 1);
//        complement(seq_forward, seq_forwardcomplement);
//
//        _gtf_feature_5prime3prime.emplace_back(seq_forward);
//        _gtf_seq_forwardcomplement.emplace_back(seq_forwardcomplement);
//    }
//
//    std::cout << "Found " << _gtf_data.size() << " features..." << std::endl;
//}

void build_index::extract_features(const std::string& feature, bool noalt) {
    // gets 5' to 3' sequence of all features

    // split gtf lines on tabs
    // gtf format
    // 0: contig, 1: source, 2: feature, 3: start, 4: end, 5: score, 6: strand, 7: frame, 8: geneID

    for (auto & l : _gtf_lines){
        std::istringstream ss(l);
        std::string substring;
        std::tuple<std::string,
                std::string,
                std::string,
                std::uint32_t,
                std::uint32_t,
                std::string,
                std::string,
                std::string,
                std::string> row_tuple =  std::make_tuple("", "", "", 0, 0, "", "", "", "");
        int i = 0;
        while (getline(ss, substring, '\t')){
            // TODO replace with case switch or just avoid this entirely
            if (i==0){
                std::get<0>(row_tuple) = substring;
            } else if (i==1) {
                std::get<1>(row_tuple) = substring;
            } else if (i==2) {
                std::get<2>(row_tuple) = substring;
            } else if (i==3) {
                uint32_t x = std::stoi(substring);
                std::get<3>(row_tuple) = x;
            } else if (i==4) {
                uint32_t x = std::stoi(substring);
                std::get<4>(row_tuple) = x;
            } else if (i==5) {
                std::get<5>(row_tuple) = substring;
            } else if (i==6) {
                std::get<6>(row_tuple) = substring;
            } else if (i==7) {
                std::get<7>(row_tuple) = substring;
            } else if (i==8) {
                std::get<8>(row_tuple) = substring;
            }
            i++;
        }

        // skip alt scaffolds from gtf if requested
        if (noalt){
            std::string contig;
            bool alt = false;
            contig = std::get<0>(row_tuple);
            if (contig.find("alt") != std::string::npos){
                alt = true;
            } //else if (contig.find('_') != std::string::npos) {
              //  alt = true;
            //}
            if (alt){
                continue;
            }
        }

        // only keep rows with feature type matching given feature
        if (std::get<2>(row_tuple) == feature){
            _gtf_data.emplace_back(row_tuple);
        }
    }

    // get sequence on gene strand (forward or reverse depending on gene) 5' to 3' for each gene
    std::string contigname;
    size_t contigname_idx;
    uint32_t start;
    uint32_t end;
    std::string strand;
    std::string seq_forward;
    std::string seq_forwardcomplement;
    std::string geneID;
    for (auto & x: _gtf_data){
        contigname = std::get<0>(x);
        start = std::get<3>(x);
        end = std::get<4>(x);
        strand = std::get<6>(x);

        auto it = std::find(_contigname_vec.begin(), _contigname_vec.end(), contigname);
        if (it == _contigname_vec.end()){
            std::cout << "contig " << contigname << " from annotation not found in reference fasta, consider removing alt scaffolds with --no-alt" << std::endl;
            continue;
        }

        contigname_idx = it - _contigname_vec.begin();
        seq_forward = _seq_vec[contigname_idx].substr(start, end - start + 1);

        if (strand == "+"){
            _gtf_feature_5prime3prime.emplace_back(seq_forward);
        } else if (strand == "-"){
            // complement and reverse to get reverse complement
            complement(seq_forward, seq_forwardcomplement);
            std::reverse(seq_forwardcomplement.begin(), seq_forwardcomplement.end());
            _gtf_feature_5prime3prime.emplace_back(seq_forwardcomplement);
        } else{
            std::cout << "strand not + or - in gtf file" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    std::cout << "Found " << _gtf_data.size() << " features..." << std::endl;
}

void build_index::build(int k, int m) {
    // build new index
    // TODO multithreaded building is possible

    // get k
    _k = k;

    // _table maps nucleotide kmer to 64 bit int containing all needed information
    std::string kmer;
    std::tuple<std::string, std::uint32_t, bool, std::string> kmer_info;
    std::vector<uint64_t> * kmer_hit_info_vec;
    std::string contig;
    uint32_t seq_start_position;
    uint32_t seq_end_position;
    std::string strand;
    std::string geneID;


    for (uint32_t i=0; i < _gtf_data.size(); ++i){
        std::cout << i << " of " << _gtf_data.size() << std::endl;
        // insert kmer 5' to 3' on forward strand
        // position is of 5' end
        std::string *seq_feature = &_gtf_feature_5prime3prime[i];
        contig = std::get<0>(_gtf_data[i]);
        seq_start_position = std::get<3>(_gtf_data[i]);
        seq_end_position = std::get<4>(_gtf_data[i]);
        strand = std::get<6>(_gtf_data[i]);
        geneID = std::get<8>(_gtf_data[i]);

        // TODO why is anything ever generating null sequence? single nucleotide exon?

        // avoid crashing on bad seq
        if (*(void**)seq_feature != nullptr){ // check NULL
            if (seq_feature->size() <= k){ // skip sequences below length of kmer
                continue;
            }
        } else {
            continue;
        }

        for (uint32_t j=0; j < seq_feature->size() - k; ++j){
            kmer = seq_feature->substr(j, k);
            // create 64 bit int representing info
            uint64_t kmer_bitset = 0;
            // 0:16 : geneID
            uint64_t geneIDidx;
            auto geneIDpos = std::find(bm.geneID_array.begin(), bm.geneID_array.end(), geneID);
            if (geneIDpos != bm.geneID_array.end()){
                // use geneID index if already added
                geneIDidx = geneIDpos - bm.geneID_array.begin();
            } else {
                // add new geneID
                bm.geneID_array.emplace_back(geneID);
                geneIDidx = bm.geneID_array.size() - 1;
            }
            bm.setbits(kmer_bitset, geneIDidx, 0, 16);
            // 16 : strand, stays set to 0 if on reverse strand
            if (strand == "+"){
                bm.setbits(kmer_bitset, (uint64_t)1, 16, 17);
            } else if (strand != "-"){
                std::cout << "strand not + or - in gtf file" << std::endl;
                exit(EXIT_FAILURE);
            }
            // 17:32 : contig
            uint64_t contigidx;
            auto contigpos = std::find(bm.contig_array.begin(), bm.contig_array.end(), contig);
            if (contigpos != bm.contig_array.end()){
                // use contig index if already added
                contigidx = contigpos - bm.contig_array.begin();
            } else{
                // add new geneID
                bm.contig_array.emplace_back(contig);
                contigidx = bm.contig_array.size() - 1; //TODO ensure not obo
            }
            bm.setbits(kmer_bitset, contigidx, 17, 32);
            // 32:64 : position within contig
            // CAREFUL think about genes on reverse complement
            uint64_t pos;
            if (strand == "+"){
                pos = seq_start_position + j; // position of first 5' nucleotide in kmer
                bm.setbits(kmer_bitset, pos, 32, 64);
            } else if (strand == "-"){
                pos = seq_end_position - j; // position of first 5' nucleotide in kmer
                bm.setbits(kmer_bitset, pos, 32, 64);
            } else{
                std::cout << "strand not + or - in gtf file" << std::endl;
                exit(EXIT_FAILURE);
            }

            // add bitset to table
            kmer_hit_info_vec = & _table[kmer];
            if (kmer_hit_info_vec->size() < m){
                kmer_hit_info_vec->emplace_back(kmer_bitset);
            }
        }
    }
}

//void showtime(const char *name, std::function<void ()> doit)
//{
//    auto t1 = std::chrono::high_resolution_clock::now();
//    doit();
//    auto t2 = std::chrono::high_resolution_clock::now();
//    auto elapsed = milliseconds<double>(t2 - t1).count();
//    printf("%s: %.3fs\n", name, (int)elapsed / 1000.0f);
//}

void build_index::serialize_table(const std::string& path_out) {
    std::cout << "Saving index to " << path_out << std::endl;

    size_t s = _table.size();
    std::cout << s << " unique kmers" << std::endl;

//    showtime("build hash", [this, path_out]() {
    std::ofstream os(path_out, std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
//}
    archive(_k);
    archive(s);
    archive(_table);

    archive(bm.geneID_array);
    archive(bm.contig_array);
}


void build_index::complement(std::string &s, std::string &comp) {
    std::string c;
    for (char ch : s){
        if (ch == 'A'){
            c += 'T';
        } else if (ch == 'T'){
            c += 'A';
        } else if (ch == 'G'){
            c += 'C';
        } else if (ch == 'C'){
            c += 'G';
        } else{
            c += ch; // just preserve non-ATGC nucleotides for now
        }
        // TODO nonstandard nucleotide table
    }
    comp = c;
}




