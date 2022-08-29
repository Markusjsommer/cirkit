#include <iostream>
#include "cxxopts.hpp"
#include "FastaReader.h"
#include "build_index.h"
#include <iterator>

//#include <parallel_hashmap/phmap.h>
#include <kmer_flathash.h>
#include <bitset>
#include "BitsetManager.h"

struct compare_genevec {
    // comparison of the first value of two gene tuples
    // used to check if kmers are in the same gene locus
    bool operator()(std::tuple<uint32_t, bool, uint32_t> a, std::tuple<uint32_t, bool, uint32_t> b) const {
        return std::get<0>(a) == std::get<0>(b);
    };
};

void complement(std::string &s, std::string &comp) {
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


int main(int argc, char* argv[]) {
    // parse command line options
    cxxopts::Options options("cirkit", "Fast, conservative, reference-guided circRNA identification");
    options.add_options()
            ("1, reads-1", "Path to input fasta for paired end mate 1", cxxopts::value<std::string>())
            ("2, reads-2", "Path to input fasta for paired end mate 2", cxxopts::value<std::string>())
            ("o, out", "Path to directory to output circRNA info", cxxopts::value<std::string>())
            ("index", "Path of kmer index", cxxopts::value<std::string>())

            ("build", "Build index from GTF and reference fasta, requires --gtf, --ref, and --index", cxxopts::value<bool>()->default_value("false"))
            ("gtf", "Path to GTF with transcript annotation", cxxopts::value<std::string>())
            ("ref", "Path to reference fasta for building index", cxxopts::value<std::string>())
            ("feature", "Feature name in GTF where circRNAs should be found, e.g. CDS, locus, transcript, exon; preprocessing an annotation with gffread -M may help discard intergenic regions and redundant transcripts", cxxopts::value<std::string>()->default_value("CDS"))
            ("no-alt", R"(Filter out alternate scaffolds from reference (really just removes any contigs with "alt" or "fix" in the name))", cxxopts::value<bool>()->default_value("false"))
            ("k", "nucleotide kmer length", cxxopts::value<int>()->default_value("25"))

//            ("k", "nucleotide kmer length", cxxopts::value<int>()->default_value("35")) // TODO minimizers
//            ("l", "minimizer length", cxxopts::value<int>()->default_value("31"))
            ("m, max-mappers", "maximum number of kmer locations to save in index, lower number can stop degenerate sequence from being too much of a problem", cxxopts::value<int>()->default_value("10"))
            ("h, help", "Print cirkit usage")
            ;

    auto result = options.parse(argc, argv);

    // check validity and display help
    if (result.count("help") or (not (result.count("build") or (result.count("1") and result.count("2") and result.count("out") and result.count("index"))))){
        std::cout << options.help() << std::endl;
        return 0;
    }

    // build index if requested
    if (result.count("build")){
        std::cout << "Building index..." << std::endl;

        // check for required args
        if (not (result.count("gtf") and (result.count("ref")) and (result.count("index")))){
            std::cout << "Must provide --gtf, --ref, and --index if building index" << std::endl;
            return 1;
        }

        // read reference and build new index
        bool noalt = result["no-alt"].as<bool>();
        std::string gtf = result["gtf"].as<std::string>();
        std::string ref = result["ref"].as<std::string>();
        std::string table_out = result["index"].as<std::string>();
        std::string feature = result["feature"].as<std::string>();
        build_index index;
        index.read_ref(gtf, ref);
        index.extract_features(feature, noalt);
        index.build(result["k"].as<int>(), result["m"].as<int>());
        index.serialize_table(table_out);

        std::cout << "Done" << std::endl;
        return 0;
    }

    // load input fastas
    // TODO replace with a real fasta/fastq reader
    FastaReader fastaio;
    std::vector<std::string> seq_vec_r1;
    std::vector<std::string> contigname_vec_r1;
    std::cout << "Reading from " << result["1"].as<std::string>() << std::endl;
    int fasta_fail = fastaio.read_fasta(result["1"].as<std::string>(), seq_vec_r1, contigname_vec_r1);
    if (fasta_fail) {
        return 1;
    }
    std::vector<std::string> seq_vec_r2;
    std::vector<std::string> contigname_vec_r2;
    std::cout << "Reading from " << result["2"].as<std::string>() << std::endl;
    fasta_fail = fastaio.read_fasta(result["2"].as<std::string>(), seq_vec_r2, contigname_vec_r2);
    if (fasta_fail){
        return 1;
    }

    // check number of paired reads
    if (seq_vec_r1.size() != seq_vec_r2.size()){
        std::cout << "Paired read files contain different numbers of reads" << std::endl;
        return 1;
    }

    // capitalize all reads
    for (auto &seq : seq_vec_r1){
        for(auto &c: seq){
            c = toupper(c);
        }
    }
    for (auto &seq : seq_vec_r2){
        for(auto &c: seq){
            c = toupper(c);
        }
    }

    // load index
    std::cout << "Loading index from " << result["index"].as<std::string>() << std::endl;
    kmer_flathash table;
    table.load(result["index"].as<std::string>());

    // get size of k used for index
    int k = table._k;
    std::cout << "Index uses kmer length "  << k << std::endl;






    // look for circRNAs
    std::cout << "Checking for circRNA kmer signatures..." << std::endl;
    std::vector<std::vector<uint64_t>> putative_circRNA_kmers;
    std::vector<size_t> putative_circRNA_readidx;

    std::string *seq_m1;
    std::string *seq_m2;

    std::string kmer_m1;
    std::string kmer_m2;
    std::string kmer_m1_rc;
    std::string kmer_m2_rc;

    BitsetManager bm;

    size_t n_mapped_linear = 0;
    size_t n_mapped_circular = 0;
    size_t n_unmapped = 0;
    size_t n_mapped_exactloc_FOOOO = 0; //TODO is this more than expected?
    std::vector<std::tuple<std::string, size_t, size_t, bool>> pair_map_info_circ_vec; // for mapping read pair <contig, coord r1, coord r2, strand>
    std::vector<std::tuple<std::string, size_t, size_t, bool>> pair_map_info_linear_vec; // for mapping read pair <contig, coord r1, coord r2, strand>

    for (size_t i=0; i < seq_vec_r1.size(); ++i){
//        std::cout << i+1 << " of "  << seq_vec_r1.size() << std::endl;

        seq_m1 = &seq_vec_r1[i];
        seq_m2 = &seq_vec_r2[i];

        // start from 5' kmer on mate 1 and 5' kmer on mate 2 (opposite ends of read pair)


        // search inwards from 5' ends of each pair for first match
        // check concordant at same locus
            // break
        // check discordant at same locus
            // add to putative circRNA list

        // how to filter for good circRNAs?
        // if 1 read is split mapped (read crosses bsj) then second read should map to same locus
        // don't use unplaced contigs
        // mask repetitive regions (still a problem in exons only?)
        // more than 1 read spanning same bsj
        // insert size of circRNA should appear relatively large



        bool decided = false;
        bool mapped = false;

        std::vector<std::tuple<uint32_t, bool, uint32_t>> genevec_mapping_m1; // store all (geneint, strand, position) for each read to check for shared mappings
        std::vector<std::tuple<uint32_t, bool, uint32_t>> genevec_mapping_m2;
        std::vector<std::tuple<uint32_t, bool, uint32_t>> genevec_mapping_m1_rc;
        std::vector<std::tuple<uint32_t, bool, uint32_t>> genevec_mapping_m2_rc;

        size_t kmer_idx_m1 = 0;
        size_t kmer_idx_m2 = 0;
        bool in_table;
        bool in_table_m1;
        bool in_table_m1_rc;
        bool in_table_m2;
        bool in_table_m2_rc;

        bool mapped_strand_m1;
        bool mapped_strand_m2;
        size_t mapped_loc_m;
        size_t mapped_loc_m_rc;
        std::string mapped_name;

        std::vector<uint64_t> table_value;
        while (not decided){
            // find farthest kmer map in read pair that maps to same locus
            // search for map in read 1, from 5' end
            while (kmer_idx_m1 < seq_m1->length() - k){
                // check forward strand
                kmer_m1 = seq_m1->substr(kmer_idx_m1, k);
                in_table_m1 = table.query(kmer_m1, table_value);
                if (in_table_m1){
                    in_table = true;
                    // add geneIDs, strand, and position to running set for R1
                    for (uint64_t x : table_value) {
                        uint32_t geneint = bm.getbits(x, 0, 16);
                        bool strand = bm.getbits(x, 16, 17);
                        uint32_t pos = bm.getbits(x, 32, 64);
                        std::tuple<uint32_t, bool, uint32_t> info = std::make_tuple(geneint, strand, pos);
                        genevec_mapping_m1.emplace_back(info);
                    }
                }
                // check reverse complement strand
                complement(kmer_m1, kmer_m1_rc);
                std::reverse(kmer_m1_rc.begin(), kmer_m1_rc.end());
                in_table_m1_rc = table.query(kmer_m1_rc, table_value);
                if (in_table_m1_rc){
                    in_table = true;
                    // add geneIDs, strand, and position to running set for R1
                    for (uint64_t x : table_value) {
                        uint32_t geneint = bm.getbits(x, 0, 16);
                        bool strand = bm.getbits(x, 16, 17);
                        uint32_t pos = bm.getbits(x, 32, 64);
                        std::tuple<uint32_t, bool, uint32_t> info = std::make_tuple(geneint, strand, pos);
                        genevec_mapping_m1_rc.emplace_back(info);
                    }
                }

                // update kmer idx and proceed to checking for matches in mate 2
                kmer_idx_m1 ++;
                if (in_table){
                    in_table = false;
                    break;
                }
            }

            // search for map in read 2
            while (kmer_idx_m2 < seq_m2->length() - k){
                // check forward strand
                kmer_m2 = seq_m2->substr(kmer_idx_m2, k);
                in_table_m2 = table.query(kmer_m2, table_value);
                if (in_table_m2){
                    in_table = true;
                    // add geneIDs to running set for R2
                    for (uint64_t x : table_value) {
                        uint32_t geneint = bm.getbits(x, 0, 16);
                        bool strand = bm.getbits(x, 16, 17);
                        uint32_t pos = bm.getbits(x, 32, 64);
                        std::tuple<uint32_t, bool, uint32_t> info = std::make_tuple(geneint, strand, pos);
                        genevec_mapping_m2.emplace_back(info);
                    }
                }

                // check reverse complement strand
                complement(kmer_m2, kmer_m2_rc);
                std::reverse(kmer_m2_rc.begin(), kmer_m2_rc.end());
                in_table_m2_rc = table.query(kmer_m2_rc, table_value);
                if (in_table_m2_rc){
                    in_table = true;
                    // add geneIDs, strand, and position to running set for R1
                    for (uint64_t x : table_value) {
                        uint32_t geneint = bm.getbits(x, 0, 16);
                        bool strand = bm.getbits(x, 16, 17);
                        uint32_t pos = bm.getbits(x, 32, 64);
                        std::tuple<uint32_t, bool, uint32_t> info = std::make_tuple(geneint, strand, pos);
                        genevec_mapping_m2_rc.emplace_back(info);
                    }
                }

                // update kmer idx and proceed to checking for hits
                kmer_idx_m2 ++;
                if (in_table){
                    in_table = false;
                    break;
                }

            }

            // check if any kmer mappings overlap on geneID
            // m1 with m2_rc
            if ((not genevec_mapping_m1.empty()) and (not genevec_mapping_m2_rc.empty())){
                auto it_m1 = std::find_first_of (genevec_mapping_m1.begin(), genevec_mapping_m1.end(),
                                                 genevec_mapping_m2_rc.begin(), genevec_mapping_m2_rc.end(),
                                                 compare_genevec{});
                if (it_m1 != genevec_mapping_m1.end()){


                    auto it_m2_rc = std::find_if(genevec_mapping_m2_rc.begin(),  genevec_mapping_m2_rc.end(), [it_m1](std::tuple<uint32_t, bool, uint32_t> const& tup){ return std::get<0>(tup) == std::get<0>(*it_m1); });

//                    std::cout << table.geneID_array[std::get<0>(*it_m1)] << std::endl;
//                    std::cout << std::get<1>(*it_m1) << std::endl;
//                    std::cout << std::get<2>(*it_m1) << std::endl;
//                    std::cout << table.geneID_array[std::get<0>(*it_m2_rc)] << std::endl;
//                    std::cout << std::get<1>(*it_m2_rc) << std::endl;
//                    std::cout << std::get<2>(*it_m2_rc) << std::endl;
//                    std::cout << std::endl;


                    mapped_strand_m1 = std::get<1>(*it_m1);
                    mapped_strand_m2 = std::get<1>(*it_m2_rc);
                    // ensure both hits map to same stranded feature
                    if (mapped_strand_m1 == mapped_strand_m2){
                        mapped = true;
                        decided = true;
                        mapped_loc_m = std::get<2>(*it_m1);
                        mapped_loc_m_rc = std::get<2>(*it_m2_rc);
                        mapped_name = table.geneID_array[std::get<0>(*it_m1)];
                        continue; // skip final position check to avoid false-false-mappings
                    }
                }
            }

            // m2 and m1_rc
            if ((not genevec_mapping_m1_rc.empty()) and (not genevec_mapping_m2.empty())){
                auto it_m1_rc = std::find_first_of (genevec_mapping_m1_rc.begin(), genevec_mapping_m1_rc.end(),
                                                 genevec_mapping_m2.begin(), genevec_mapping_m2.end(),
                                                 compare_genevec{});
                if (it_m1_rc != genevec_mapping_m1_rc.end()){


                    auto it_m2 = std::find_if(genevec_mapping_m2.begin(),  genevec_mapping_m2.end(), [it_m1_rc](std::tuple<uint32_t, bool, uint32_t> const& tup){ return std::get<0>(tup) == std::get<0>(*it_m1_rc); });


//                    std::cout << table.geneID_array[std::get<0>(*it_m1_rc)] << std::endl;
//                    std::cout << std::get<1>(*it_m1_rc) << std::endl;
//                    std::cout << std::get<2>(*it_m1_rc) << std::endl;
//                    std::cout << table.geneID_array[std::get<0>(*it_m2)] << std::endl;
//                    std::cout << std::get<1>(*it_m2) << std::endl;
//                    std::cout << std::get<2>(*it_m2) << std::endl;
//                    std::cout << std::endl;

                    mapped_strand_m1 = std::get<1>(*it_m1_rc);
                    mapped_strand_m2 = std::get<1>(*it_m2);
                    // ensure both hits map to same stranded feature
                    if (mapped_strand_m1 == mapped_strand_m2){
                        mapped = true;
                        decided = true;
                        mapped_loc_m = std::get<2>(*it_m2);
                        mapped_loc_m_rc = std::get<2>(*it_m1_rc);
                        mapped_name = table.geneID_array[std::get<0>(*it_m1_rc)]; // both necessarily map to same name
                        continue; // skip final position check to avoid false-false-mappings
                    }
                }
            }

            if (kmer_idx_m1 == seq_m1->length()-k and kmer_idx_m2 == seq_m2->length()-k){
                mapped = false;
                decided = true;
            }

        }

        // distinguish linear rna from putative circRNA
        // linear rna mapping:    5' ----->          3'
        //                        3'          <----- 5'
        // circ rna mapping:      5' <-----          3'
        //                        3'          -----> 5'
        // i.e. given index of 5' to 3' gene sequence
        // linear has the reverse complement match closer to 3'
        // circ has the reverse complement match closer to 5'
        if (mapped) {
            std::tuple<std::string, size_t, size_t, bool> pair_map_info;
            if (mapped_strand_m1){
                // check for linear rna
                if (mapped_loc_m < mapped_loc_m_rc){
                    ++ n_mapped_linear;
                    pair_map_info = std::make_tuple(mapped_name,
                                                    mapped_loc_m,
                                                    mapped_loc_m_rc,
                                                    mapped_strand_m1);
                    pair_map_info_linear_vec.emplace_back(pair_map_info);
                // check for putative circular RNA
                } else if (mapped_loc_m > mapped_loc_m_rc){
                    ++ n_mapped_circular;
                    pair_map_info = std::make_tuple(mapped_name,
                                                    mapped_loc_m,
                                                    mapped_loc_m_rc,
                                                    mapped_strand_m1);
                    pair_map_info_circ_vec.emplace_back(pair_map_info);
                } else{
                    ++n_mapped_exactloc_FOOOO; // TODO remove
                }
            } else{
                // check for linear rna
                if (mapped_loc_m > mapped_loc_m_rc){
                    ++ n_mapped_linear;
                    pair_map_info = std::make_tuple(mapped_name,
                                                    mapped_loc_m,
                                                    mapped_loc_m_rc,
                                                    mapped_strand_m1);
                    pair_map_info_linear_vec.emplace_back(pair_map_info);
                // check for putative circular RNA
                } else if (mapped_loc_m < mapped_loc_m_rc){
                    ++ n_mapped_circular;
                    pair_map_info = std::make_tuple(mapped_name,
                                                    mapped_loc_m,
                                                    mapped_loc_m_rc,
                                                    mapped_strand_m1);
                    pair_map_info_circ_vec.emplace_back(pair_map_info);
                } else{
                    ++n_mapped_exactloc_FOOOO; // TODO remove
                }
            }
        } else {
            ++ n_unmapped;
        }

    }

    std::cout << "n_mapped_linear: " << n_mapped_linear << std::endl;
    std::cout << "n_mapped_circular: " << n_mapped_circular << std::endl;
    std::cout << "n_mapped_exactloc_FOOOO: " << n_mapped_exactloc_FOOOO << std::endl;
    std::cout << "n_unmapped: " << n_unmapped << std::endl;

    std::cout << n_mapped_linear + n_mapped_circular + n_unmapped + n_mapped_exactloc_FOOOO << std::endl;
    std::cout << seq_vec_r1.size() << std::endl;

    // write output files
    // linear
    std::string p_out_linear = result["o"].as<std::string>();;
    if (p_out_linear.back() == '/'){ // remove trailing slash if present
        p_out_linear.pop_back();
    }
    p_out_linear = p_out_linear + "/read_map_linearRNA.csv";
    std::cout << "Writing to " << p_out_linear << std::endl;
    std::ofstream out_linear;
    out_linear.open(p_out_linear);
    out_linear << "locus" << ","
               << "coord_m1" << ","
               << "coorc_m2" << ","
               << "strand" << "\n";
    for (auto x : pair_map_info_linear_vec){
        out_linear << std::get<0>(x) << ","
                   << std::get<1>(x) << ","
                   << std::get<2>(x) << ","
                   << std::get<3>(x) << "\n";
    }
    out_linear.close();
    // circRNA
    std::string p_out_circ = result["o"].as<std::string>();
    if (p_out_circ.back() == '/'){ // remove trailing slash if present
        p_out_circ.pop_back();
    }
    p_out_circ = p_out_circ + "/read_map_circRNA.csv";
    std::cout << "Writing to " << p_out_circ << std::endl;
    std::ofstream out_circ;
    out_circ.open(p_out_circ);
    out_circ << "locus" << ","
             << "coord_m1" << ","
             << "coorc_m2" << ","
             << "strand" << "\n";
    for (auto x : pair_map_info_circ_vec){
        out_circ << std::get<0>(x) << ","
                 << std::get<1>(x) << ","
                 << std::get<2>(x) << ","
                 << std::get<3>(x) << "\n";
    }
    out_circ.close();


    // TODO write read IDs as well, also contig ID
    // TODO figure out how to report actual bsj location rather than read mapping location
    // TODO perhaps it is best to do that in a separate step by looking at exon junctions

    // TODO filter out false positives from mapping linear genome to index

    // TODO filter out ribosomal RNA mapping

    // TODO check if out directory exists and make it if not




    return 0;
}
