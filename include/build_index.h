//
// Created by marku on 10/27/2021.
//
#include <iostream>
#include "FastaReader.h"
#include "FileReader.h"
//#include <parallel_hashmap/phmap.h>
#include <set>
#include <unordered_map>
#include <bitset>
#include "BitsetManager.h"

//using phmap::flat_hash_map;

#ifndef NEWCIRC_INDEXER_H
#define NEWCIRC_INDEXER_H


class build_index {
public:
    build_index() = default;
    ~build_index() = default;

    void read_ref(std::string &gtf, std::string &ref);
    void extract_features(const std::string& feature, bool noalt);
    void build(int k, int m);
    void serialize_table(const std::string &path_out);

private:
    static void complement(std::string &s, std::string &comp);
    FastaReader _io;
    std::vector<std::string> _seq_vec;
    std::vector<std::string> _contigname_vec;
    int _k;

    FileReader _gtfio;
    std::vector<std::string> _gtf_lines;

    // gtf format (WITH geneID replacing attributes at end to make it easier to track)
    // 0: contig, 1: source, 2: feature, 3: start, 4: end, 5: score, 6: strand, 7: frame, 8: geneID
    std::vector<std::tuple<std::string,
                           std::string,
                           std::string,
                           std::uint32_t,
                           std::uint32_t,
                           std::string,
                           std::string,
                           std::string,
                           std::string>> _gtf_data;

    // nucleotide sequence from both strands
    std::vector<std::string> _gtf_feature_5prime3prime;

    // hash table
    using MapType = std::unordered_map<std::string, std::vector<uint64_t>>;
    MapType _table;

private:
    BitsetManager bm;

};


#endif //NEWCIRC_INDEXER_H
