//
// Created by markus on 11/18/2021.
//

#include <iostream>
#include "BitsetManager.h"
#include <bitset>


void BitsetManager::setbits(uint64_t &bits, const uint64_t value, const uint_fast8_t startbit, const uint_fast8_t endbit) {
    // pythonish indexing, 0 indexing, includes startbit and goes up to but not including endbit
    uint64_t rightbits = bits & (~0ULL >> (endbit)); // zero everything left of endbit
    uint64_t leftbits = bits & (~0ULL << (64 - startbit)); // zero everything right of startbit
    uint64_t valbits = value << (64 - (endbit - startbit)); // move value to correct position and zero everything left of startbit and right of endbit
    valbits = valbits >> startbit;
    bits = leftbits | valbits | rightbits; // bitwise or to set bits
}


uint64_t BitsetManager::getbits(const uint64_t bits, const uint_fast8_t startbit, const uint_fast8_t endbit) {
    // pythonish indexing, 0 indexing, includes startbit and goes up to but not including endbit
    uint64_t value = bits;
    value = value >> (64 - endbit);
    value = value & (~0ULL >> (64 - (endbit - startbit)));
    return value;
}


//uint64_t BitsetManager::tuple_to_uint64_t(std::tuple<std::string, std::uint32_t, bool, std::string> x) {
//    // tuple is <contig, position, strand, geneID>
//    // 64 bitset for each kmer, python indexing style
//    // 0:16 : geneID
//    // 16 : strand
//    // 17:32 : contig
//    // 32:64 : position within contig
//    uint64_t bits;
//
//    return bits;
//}

//std::tuple<std::string, std::uint32_t, bool, std::string> BitsetManager::uint64_t_to_tuple(uint64_t x) {
//    // tuple is <contig, position, strand, geneID>
//
//    // 64 bitset for each kmer, python indexing style
//    // 0:16 : geneID
//    // 16 : strand
//    // 17:32 : contig
//    // 32:64 : position within contig
//    std::tuple<std::string, std::uint32_t, bool, std::string> info;
//
//
//    return info;
//}


