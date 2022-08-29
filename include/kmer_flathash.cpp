//
// Created by markus on 11/5/2021.
//


#include "kmer_flathash.h"


void kmer_flathash::load(const std::string& path) {
    std::ifstream is(path, std::ios::binary);
    cereal::BinaryInputArchive archive_in(is);
    size_t table_size;

    archive_in(_k);
    archive_in(table_size);
    _table.reserve(table_size);
    archive_in(_table);

    archive_in(geneID_array);
    archive_in(contig_array);

    archive_in(junction_info);
    for(auto kv : junction_info) {
        std::cout << kv.first << std::endl;
        std::cout << kv.second.size() << std::endl;
    }

}

bool kmer_flathash::query(const std::string& kmer, std::vector<uint64_t> &table_value) {
    auto search = _table.find(kmer);
    if (search != _table.end()){
        table_value = search->second;
        return true;
    } else{
        return false;
    }
}



