//
// Created by markus on 11/5/2021.
//

#ifndef NEWCIRC_KMER_FLATHASH_H
#define NEWCIRC_KMER_FLATHASH_H


#include <vector>
#include <iostream>
#include <fstream>
#include <iostream>
#include <bitset>
#include <cinttypes>
//#include "parallel_hashmap/phmap_dump.h"
//#include <parallel_hashmap/phmap.h>



#include <cereal/archives/binary.hpp>
#include "cereal/types/unordered_map.hpp" // this seems to cause the weird error
#include "cereal/types/memory.hpp"
#include "cereal/types/bitset.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/tuple.hpp"

#include "cereal/types/utility.hpp" // provides serialization of std::pair

#include <unordered_map>

//using phmap::flat_hash_map;


class kmer_flathash {
public:
    kmer_flathash() = default;
    ~kmer_flathash() = default;

    void load(const std::string& path);

    bool query(const std::string& kmer, std::vector<uint64_t> &table_value);
    int _k;

    std::vector<std::string> contig_array;
    std::vector<std::string> geneID_array;

    // intron-exon junction information (enables backsplice junction annotation)
    // for each contig, stores a vector of start/end pairs from the annotated exons
    std::unordered_map<std::string, // contig
            std::vector<std::pair<std::uint32_t, // start coord
                    std::uint32_t>>> junction_info; // end coord

private:
    using MapType = std::unordered_map<std::string, std::vector<uint64_t>>;
    MapType _table;



};


#endif //NEWCIRC_KMER_FLATHASH_H
