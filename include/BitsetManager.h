//
// Created by markus on 11/18/2021.
//

#ifndef CIRCUIT_BITSETMANAGER_H
#define CIRCUIT_BITSETMANAGER_H

#include <bitset>
#include <tuple>
#include <vector>

class BitsetManager {
public:
    BitsetManager() = default;
    ~BitsetManager() = default;

    // set bits
    void setbits(uint64_t &bits, uint64_t value, uint_fast8_t startbit, uint_fast8_t endbit);

    // read bits
    uint64_t getbits(uint64_t bits, uint_fast8_t startbit, uint_fast8_t endbit);

    // tuple is <contig, position, strand, geneID>

    // 64 bitset for each kmer, python indexing style
    // 0:16 : geneID
    // 16 : strand
    // 17:32 : contig
    // 32:64 : position within contig

//    uint64_t tuple_to_uint64_t(std::tuple<std::string, std::uint32_t, bool, std::string> x);
//    std::tuple<std::string, std::uint32_t, bool, std::string> uint64_t_to_tuple(uint64_t x);

    // bits to contig and geneIDs
    // position in vector is the int of related bits
    std::vector<std::string> contig_array;
    std::vector<std::string> geneID_array;

    // TODO load contig and geneID array into index serialization


};


#endif //CIRCUIT_BITSETMANAGER_H
