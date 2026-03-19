//
// Created by Patrick Flege on 16/03/2026.
//

#pragma once
#include "cli/GraphWrapper.hpp"
#include "external/KmerType.h"

#include <regex>

using node_index_kmer = GraphWrapper::node_index_kmer;
using kmer_frequencies = GraphWrapper::kmer_frequencies;
class KmerClassifier {
public:
    typedef std::vector<std::pair<int, int>> genomes_and_frequencies;
    KmerClassifier(GraphWrapper& graph_ptr);
    KmerClassifier(GraphWrapper& graph_ptr, int core_threshold, int unique_threshold);
    unsigned long long get_node_number();
    kmer_frequencies get_kmer_frequencies(std::vector<node_index_kmer>& nodes);
    int get_num_genomes();

    int extractGenomeNumber(std::string genome_contig_path);

    template<size_t genome_number> std::array<int, genome_number> count_kmer_per_genome(node_index_kmer idx) {
        std::array<int, genome_number> kmer_counts({0});
        //std::vector<int> kmer_counts(genome_number, 0);
        auto nodes = std::vector<node_index_kmer>({idx});
        auto frequencies = this->get_kmer_frequencies(nodes);
        for (size_t i = 0; i < frequencies.size(); i++) {
            auto genome_name =  std::get<0>(frequencies[i]);
            auto genome_no = this->extractGenomeNumber(genome_name);
            if (genome_no == -1)
                throw std::invalid_argument("Invalid genome number");
            auto genome_count = std::get<2>(frequencies[i])[0];
            int genome_count_int = (int)genome_count;
            kmer_counts[genome_no] = genome_count_int;
        }
        return kmer_counts;
    };


    //std::vector<int> kmer_counts
    template<size_t genome_number>
    std::tuple<KmerType, genomes_and_frequencies> classify_kmers(std::array<int, genome_number> kmer_counts) {
        KmerType type;
        genomes_and_frequencies genome_frequencies({});
        int number_of_occuring_genomes = 0;
        for (size_t i = 0; i < kmer_counts.size(); i++) {
            if (kmer_counts[i] > 0){
                number_of_occuring_genomes++;
                auto occurrence_pair = std::make_pair(i, kmer_counts[i]);
                genome_frequencies.push_back(occurrence_pair);
            }
        }

        if (number_of_occuring_genomes == this->unique_threshold) {
            std::cout << "Working" << std::endl;
            type = KmerType::UNIQUE;
        }
        else if (number_of_occuring_genomes == this->core_threshold) {
            std::cout << "Not unique" << std::endl;
            type = KmerType::CORE;
        }
        else {
            type = KmerType::ACCESSORY;
        }

        return std::tuple<KmerType, genomes_and_frequencies>(type, genome_frequencies);
    };

    template <size_t genome_number>
    void fill_kmer_matrix(std::array<std::array<int, genome_number>, 3>& kmer_matrix,
                          KmerType type,
                          genomes_and_frequencies genome_frequencies) {
        switch (type) {
            case KmerType::CORE: {
                std::cout << "Core" << std::endl;
                for_each(genome_frequencies.begin(), genome_frequencies.end(),
                         [&kmer_matrix](auto gen_idx) {
                             std::cout << "CORE, at index: " << gen_idx.first << std::endl;
                             std::cout << "CORE, adding: " << gen_idx.second << std::endl;

                             kmer_matrix[0][gen_idx.first] += gen_idx.second;
                         });
                break;
            }
            case KmerType::ACCESSORY: {
                std::cout << "Accessory" << std::endl;
                std::for_each(genome_frequencies.begin(), genome_frequencies.end(),
                              [&kmer_matrix](auto gen_idx) {
                                  std::cout << "Accessory at index: " << gen_idx.first
                                            << std::endl;
                                  std::cout << "Accessory, adding: " << gen_idx.second
                                            << std::endl;

                                  kmer_matrix[1][gen_idx.first] += gen_idx.second;
                              });
                break;
            }
            case KmerType::UNIQUE: {
                std::cout << "Unique" << std::endl;
                std::cout << "Unique, adding: " << genome_frequencies[0].second << std::endl;
                kmer_matrix[2][genome_frequencies[0].first] += genome_frequencies[0].second;
                break;
            }
            default:
                std::cout << "Unknown" << std::endl;
                break;
        }
    }
    template <size_t genome_number>
    std::array<std::array<int, genome_number>, 3> create_kmer_classification_matrix() {
        // kmer_matrix = std::vector<std::vector<int>>(3, genome_number);
        auto kmer_matrix = std::array<std::array<int, genome_number>, 3>{{{}}};
        auto number_of_kmers = this->get_node_number();
        std::cout << "Node number:" << number_of_kmers <<  std::endl;
        node_index_kmer current = 1;
        for (size_t i = 0; i < kmer_matrix.size(); i++) {
            for (size_t j = 0; j < kmer_matrix[i].size(); j++) {
                kmer_matrix[i][j] = 0;
            }
        }

        while (current < number_of_kmers) {

            auto counts_per_genomes = this->count_kmer_per_genome<genome_number>(current);
            auto classified_kmers = this->classify_kmers(counts_per_genomes);
            //TODO: Continue here.
            //TODO: Row 0 is core, row 1 accessory, row 2 unique
            KmerType type = std::get<0>(classified_kmers);
            genomes_and_frequencies genome_frequencies = std::get<1>(classified_kmers);

            fill_kmer_matrix<genome_number>(kmer_matrix, type, genome_frequencies);

            current++;
        }
        return kmer_matrix;
    }

private:
    GraphWrapper& graph_ptr;
    int core_threshold;
    int unique_threshold;
};
