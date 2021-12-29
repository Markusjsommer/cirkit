//
// Created by marku on 10/27/2021.
//

#include "FileReader.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <gzip/decompress.hpp>
#include <gzip/utils.hpp>

void FileReader::read_gtf(const std::string &filepath_in, std::vector<std::string> &line_vec) {
    // read gtf
    std::ifstream ifile(filepath_in);
    if(ifile) {
        // read file buffer into string
        std::stringstream buffer;
        buffer << ifile.rdbuf();
        std::string filedata = buffer.str();

        // check if gzipped
        const char *compressed_pointer = filedata.data();
        std::size_t compressed_size = filedata.size();
        bool gzipped = gzip::is_compressed(compressed_pointer, compressed_size);
        if (gzipped) {
            // decompress and read all data into string
            std::string decompressed_data = gzip::decompress(compressed_pointer, compressed_size);

            // read compressed data
            std::istringstream istring(decompressed_data);
            std::string line;
            while (std::getline(istring, line)) {
//                if (line.back() == '\r') {
//                    line.pop_back();
//                }
//                line_vec.back() += line;
                line_vec.emplace_back(line);
            }
        } else {
            // read text fasta data
            std::istringstream istring(filedata);
            std::string line;
            while (std::getline(istring, line)) {
//                if (line.back() == '\r'){
//                    line.pop_back();
//                } else {
//                    line_vec.back() += line;
//                }
                line_vec.emplace_back(line);
            }
        }
    } else {
        std::cout << "Error opening gtf file" << std::endl;
    }
}
