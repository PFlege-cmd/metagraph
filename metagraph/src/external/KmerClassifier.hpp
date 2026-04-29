//
// Created by Patrick Flege on 16/03/2026.
//

#pragma once
#include "cli/GraphWrapper.hpp"
#include "external/KmerType.h"

#include <regex>
#include <vector>

using node_index_kmer = GraphWrapper::node_index_kmer;
using kmer_frequencies = GraphWrapper::kmer_frequencies;
class KmerClassifier {
public:
    typedef std::vector<std::pair<int, int>> genomes_and_frequencies;
    KmerClassifier(GraphWrapper& graph_ptr);
    KmerClassifier(GraphWrapper& graph_ptr, int core_threshold, int unique_threshold);
    KmerClassifier(GraphWrapper& graph_ptr, int core_threshold, int unique_threshold, int genome_number);
    unsigned long long get_node_number();
    kmer_frequencies get_kmer_frequencies(std::vector<node_index_kmer>& nodes);
    int get_num_genomes();
    void set_num_genomes(int genome_number);
    int extractGenomeNumber(std::string genome_contig_path);


    std::vector<int> count_kmer_per_genome(node_index_kmer idx) {
        std::vector<int> kmer_counts(this->genome_number, {0});
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

    std::tuple<KmerType, genomes_and_frequencies> classify_kmers(std::vector<int> kmer_counts) {
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
    /**
     *
     * @param genome_frequencies genome and frequencies of a given kmer
     * @param gen_idx pair<int, int> of index,frequencies
     */
    void fill_triangular_matrices(const genomes_and_frequencies& genome_frequencies, std::pair<int, int> gen_idx);


    /**
     * @param type Kmertype (Unique, Accessory or Core)
     * @param genome_frequencies Vector containing pairs genome-number -  frequency of kmer
     */
    void fill_kmer_matrix(KmerType type,
                          genomes_and_frequencies genome_frequencies) {
        bool* is_present_in_genome = new bool[genome_number];
        switch (type) {
            case KmerType::CORE: {
                std::cout << "Core" << std::endl;
                for_each(genome_frequencies.begin(), genome_frequencies.end(),
                         [this, genome_frequencies, is_present_in_genome](auto gen_idx) {
                             std::cout << "CORE, at index: " << gen_idx.first << std::endl;
                             std::cout << "CORE, adding: " << gen_idx.second << std::endl;
                             this->total_kmer_matrix[0][gen_idx.first] += gen_idx.second;
                             this->distinct_kmer_matrix[0][gen_idx.first] += 1;
                             this->fill_triangular_matrices(genome_frequencies, gen_idx);
                            is_present_in_genome[gen_idx.first] = true;
                         });
                break;
            }
            case KmerType::ACCESSORY: {
                std::cout << "Accessory" << std::endl;
                std::for_each(genome_frequencies.begin(), genome_frequencies.end(),
                              [this, genome_frequencies,is_present_in_genome](auto gen_idx) {
                                  std::cout << "Accessory at index: " << gen_idx.first
                                            << std::endl;
                                  std::cout << "Accessory, adding: " << gen_idx.second
                                            << std::endl;

                                  this->total_kmer_matrix[1][gen_idx.first] += gen_idx.second;
                                  this->distinct_kmer_matrix[1][gen_idx.first] += 1;
                                  this->fill_triangular_matrices(genome_frequencies, gen_idx);
                                  is_present_in_genome[gen_idx.first] = true;

                              });
                break;
            }
            case KmerType::UNIQUE: {
                std::cout << "Unique" << std::endl;
                std::cout << "Unique, adding: " << genome_frequencies[0].second << std::endl;
                this->total_kmer_matrix[2][genome_frequencies[0].first] += genome_frequencies[0].second;
                this->distinct_kmer_matrix[2][genome_frequencies[0].first] += 1;

                this->all_shared_matrix[genome_frequencies[0].first][genome_frequencies[0].first] += genome_frequencies[0].second;
                this->all_total_matrix[genome_frequencies[0].first][genome_frequencies[0].first] += genome_frequencies[0].second;

                this->distinct_shared[genome_frequencies[0].first][genome_frequencies[0].first]++;
                this->distinct_total[genome_frequencies[0].first][genome_frequencies[0].first]++;
                is_present_in_genome[genome_frequencies[0].first] = true;
                break;
            }
            default:
                std::cout << "Unknown" << std::endl;
                break;
        }

        // for_each(genome_frequencies.begin(), genome_frequencies.end(),[this, is_present_in_genome](auto gen_idx) {
        //     for (int i = 0; i < this->genome_number; i++) {
        //         if (!is_present_in_genome[i]){
        //             this->all_total_matrix[gen_idx.first][i] += gen_idx.second;
        //             this->all_total_matrix[i][gen_idx.first] += gen_idx.second;
        //
        //             ++this->distinct_total[gen_idx.first][i];
        //             ++this->distinct_total[i][gen_idx.first];
        //         }
        //     }
        // });
        std::vector<int> frequencies = std::vector<int>(this->genome_number);
        for_each(genome_frequencies.begin(), genome_frequencies.end(), [frequencies_ref = &frequencies](auto gen_idx) {
            (*frequencies_ref)[gen_idx.first] += gen_idx.second;
        });

        for (int i = 0; i < genome_number; i++) {
            for (int j = 0; j < genome_number; j++) {
                std::cout << "Frequencies are: " << frequencies[i] << std::endl;
                int highest = frequencies[i] > frequencies[j]? frequencies[i] : frequencies[j];
                //int lowest = frequencies[i] < frequencies[j]? frequencies[i] : frequencies[j];
                if (highest == 0)
                    continue;

                ++this->distinct_total[i][j];
                this->all_total_matrix[i][j] += highest;
            }
        }
    }

    /**
     * Flattens 2D matrix of vectors to a 1D int array, which is returned.
     * @param matrix 2D vector of vector<int>
     */
    int* flatten(const std::vector<std::vector<int>>& matrix);

    /**
    Creates kmer-matrices, both including total counts and distinct counts (adds only one instead of frequency).
    Stores those in instance variables:
    total_kmer_matrix, distinct_kmer_matrix
     */
    std::vector<std::vector<int>> create_kmer_classification_matrix() {

        this->total_kmer_matrix = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));
        //auto total_kmer_matrix = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));
        this->distinct_kmer_matrix = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));

        this->distinct_shared = std::vector<std::vector<int>>(this->genome_number, std::vector<int>(this->genome_number, 0));
        this->distinct_total = std::vector<std::vector<int>>(this->genome_number, std::vector<int>(this->genome_number, 0));
        this->all_shared_matrix = std::vector<std::vector<int>>(this->genome_number, std::vector<int>(this->genome_number, 0));
        this->all_total_matrix = std::vector<std::vector<int>>(this->genome_number, std::vector<int>(this->genome_number, 0));

        auto number_of_kmers = this->get_node_number();
        std::cout << "Node number:" << number_of_kmers <<  std::endl;
        node_index_kmer current = 1;
        for (size_t i = 0; i < total_kmer_matrix.size(); i++) {
            for (size_t j = 0; j < total_kmer_matrix[i].size(); j++) {
                total_kmer_matrix[i][j] = 0;
            }
        }

        while (current < number_of_kmers) {
            auto seq = this->graph_ptr.get_graph()->get_graph().get_node_sequence(current);
            std::cout << "Node sequence:" << seq << std::endl;

            auto counts_per_genomes = this->count_kmer_per_genome(current);
            auto classified_kmers = this->classify_kmers(counts_per_genomes);
            //TODO: Continue here.
            //TODO: Row 0 is core, row 1 accessory, row 2 unique
            KmerType type = std::get<0>(classified_kmers);
            genomes_and_frequencies genome_frequencies = std::get<1>(classified_kmers);
            //TODO: FIll distinct counts as well.
            fill_kmer_matrix(type, genome_frequencies);

            current++;
        }
        return total_kmer_matrix;
    }

    std::vector<std::vector<int>> get_distinct_kmer_matrix();

    std::vector<std::vector<int>> get_total_kmer_matrix();
    std::vector<std::vector<int>> get_all_total_matrix();
    std::vector<std::vector<int>> get_all_shared_matrix();
    std::vector<std::vector<int>> get_distinct_total_matrix();
    std::vector<std::vector<int>> get_distinct_shared_matrix();

  private:
    GraphWrapper& graph_ptr;
    std::vector<std::vector<int>> distinct_kmer_matrix;
    std::vector<std::vector<int>> total_kmer_matrix;

    std::vector<std::vector<int>> all_shared_matrix;
    std::vector<std::vector<int>> all_total_matrix;

    std::vector<std::vector<int>> distinct_shared;
    std::vector<std::vector<int>> distinct_total;

    int core_threshold;
    int unique_threshold;
    int genome_number;
};
