//
// Created by marku on 10/27/2021.
//

#include <vector>
#include <string>

#ifndef NEWCIRC_GTFREADER_H
#define NEWCIRC_GTFREADER_H


class FileReader {
public:
    FileReader() = default;
    ~FileReader() = default;
    void read_gtf(const std::string &filepath_in, std::vector<std::string> &line_vec);
};


#endif //NEWCIRC_GTFREADER_H
