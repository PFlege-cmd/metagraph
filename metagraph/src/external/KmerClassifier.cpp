//
// Created by Patrick Flege on 16/03/2026.
//

#include "KmerClassifier.hpp"

#include <regex>
#include <tuple>
#include <tuple>
#include <tuple>
#include <tuple>

KmerClassifier::KmerClassifier(GraphWrapper& graph_ptr) : graph_ptr(graph_ptr) {
}

 KmerClassifier::KmerClassifier(GraphWrapper& graph_ptr,
                               int core_threshold,
                               int unique_threshold) :graph_ptr(graph_ptr),  core_threshold(core_threshold), unique_threshold(unique_threshold) {
}

 KmerClassifier::KmerClassifier(GraphWrapper& graph_ptr,
                                int core_threshold,
                                int unique_threshold,
                                int genome_number) : graph_ptr(graph_ptr) , core_threshold(core_threshold), unique_threshold(unique_threshold) , genome_number(genome_number) {
 }

std::vector<std::vector<int>> KmerClassifier::get_distinct_kmer_matrix() {
    return this->distinct_kmer_matrix;
}

std::vector<std::vector<int>> KmerClassifier::get_total_kmer_matrix() {
    return this->total_kmer_matrix;
}

std::vector<std::vector<int>> KmerClassifier::get_all_total_matrix() {
    return this->all_total_matrix;
}

std::vector<std::vector<int>> KmerClassifier::get_all_shared_matrix() {
    return this->all_shared_matrix;
}

std::vector<std::vector<int>> KmerClassifier::get_distinct_total_matrix() {
    return this->distinct_total;
}

std::vector<std::vector<int>> KmerClassifier::get_distinct_shared_matrix() {
    return this->distinct_shared;
}

void KmerClassifier::set_num_genomes(int genome_number) {
    this->genome_number = genome_number;
}

unsigned long long KmerClassifier::get_node_number() {
    return graph_ptr.get_number_nodes();
}

kmer_frequencies KmerClassifier::get_kmer_frequencies(std::vector<node_index_kmer>& nodes) {
    return graph_ptr.get_kmer_frequencies(nodes);
}

int KmerClassifier::get_num_genomes() {
    return graph_ptr.get_num_genomes();
}

int KmerClassifier::extractGenomeNumber(std::string genome_contig_path) {
    std::regex rx("([0-9]+)_contig");
    std::smatch match;
    bool found = std::regex_search(genome_contig_path, match, rx);
    if (found) {
        for (std::size_t i = 0; i < match.size(); i++) {
            auto match_var = match[i].str();
            auto pos_underscore = match_var.find('_');
            match_var = match_var.substr(0, pos_underscore);
            return std::stoi(match_var);
        }
    }
    return -1;
}

void KmerClassifier::fill_triangular_matrices(const genomes_and_frequencies& genome_frequencies,
                                              std::pair<int, int> outer_genome_idx) {
    std::for_each(genome_frequencies.begin(), genome_frequencies.end(),[this, outer_genome_idx](auto inner_genome_idx) {
        int outer_index = outer_genome_idx.first;
        int inner_index = inner_genome_idx.first;

        int lowest  =  outer_genome_idx.second < inner_genome_idx.second?outer_genome_idx.second : inner_genome_idx.second;

        if (lowest > 0) {
            ++distinct_shared[outer_index][inner_index];
        }

        this->all_shared_matrix[outer_index][inner_index] += lowest;
        //this->all_total_matrix[outer_index][inner_index] += highest;
    });
}

void KmerClassifier::fill_triangular_matrices_parallel(const genomes_and_frequencies& genome_frequencies,
                                              std::pair<int, int> outer_genome_idx,
                                              std::vector<std::vector<int>>& distinct_shared,
                                              std::vector<std::vector<int>>& all_shared_matrix) {
    std::for_each(genome_frequencies.begin(), genome_frequencies.end(),[&](auto inner_genome_idx) {
        int outer_index = outer_genome_idx.first;
        int inner_index = inner_genome_idx.first;

        int lowest  =  outer_genome_idx.second < inner_genome_idx.second?outer_genome_idx.second : inner_genome_idx.second;
        //int highest =  outer_genome_idx.second > inner_genome_idx.second?outer_genome_idx.second : inner_genome_idx.second;

        if (lowest > 0) {
            ++distinct_shared[outer_index][inner_index];
        }
        //++this->distinct_total[outer_index][inner_index];

        all_shared_matrix[outer_index][inner_index] += lowest;
        //this->all_total_matrix[outer_index][inner_index] += highest;
    });
}

int* KmerClassifier::flatten(const std::vector<std::vector<int>>& matrix) {

    int* flattened = new int[genome_number*genome_number];

    for (int i = 0; i < genome_number; i++) {
        for (int j = 0; j < genome_number; j++) {
            flattened[i*genome_number + j] = matrix[i][j];
        }
    }
    return flattened;
}

std::pair<std::vector<std::vector<int>>, PangenomeStructureAdapter::PangenomeStructureCollector&> KmerClassifier::initialize_kmer_frequencies(PangenomeStructureAdapter::PangenomeStructureCollector& loop_collector, const std::vector<unsigned long long>& all_node_ids) {
    std::vector<std::vector<int>> freqs(all_node_ids.size(), std::vector<int>(genome_number, 0));

    // std::for_each(all_node_ids.begin(), all_node_ids.end(), [&](const auto& node_idx) {
    //     auto counts_per_genomes{this->count_kmer_per_genome(node_idx)};
    //     auto classified_kmers = this->classify_kmers(counts_per_genomes);
    //
    //     int reduced{std::reduce(counts_per_genomes.begin(), counts_per_genomes.end())};
    //
    //     if (reduced > 0) {
    //         freqs.push_back(counts_per_genomes);
    //     }
    //     //std::cout << "case:  " <<  kmer_type << std::endl;
    //     switch (std::get<0>(classified_kmers)) {
    //         case KmerType::CORE: {
    //             std::for_each(
    //                 loop_collector.core_loops.begin(),
    //                 loop_collector.core_loops.end(),
    //                 [&](auto& core_kmer_freq_combs)
    //                 {
    //                 core_kmer_freq_combs[this->genome_number - 2] += 1;
    //             });
    //             break;
    //         }
    //
    //         case KmerType::ACCESSORY: {
    //             std::for_each(
    //                 loop_collector.accessory_loops.begin(),
    //                 loop_collector.accessory_loops.end(),
    //                 [&](auto& accessory_kmer_freq_combs)
    //                 {
    //                 accessory_kmer_freq_combs[this->genome_number - 2] += 1;
    //              });
    //             break;
    //         }
    //
    //         case KmerType::UNIQUE: {
    //             std::for_each(
    //                 loop_collector.unique_loops.begin(),
    //                 loop_collector.unique_loops.end(),
    //                 [&](auto& unique_kmer_freq_combs)
    //                 {
    //                 unique_kmer_freq_combs[this->genome_number - 2] += 1;
    //              });
    //             break;
    //         }
    //         case KmerType::NOT_PRESENT: {
    //             break;
    //         }
    //     }
    // });

    #pragma omp parallel
    {
        #pragma omp for schedule(static, 3)
        for (node_index_kmer node_idx = 0; node_idx < all_node_ids.size(); node_idx++) {
            auto counts_per_genomes{this->count_kmer_per_genome(all_node_ids[node_idx])};
            auto classified_kmers = this->classify_kmers(counts_per_genomes);

            //int reduced{std::reduce(counts_per_genomes.begin(), counts_per_genomes.end())};

            //if (reduced > 0) {
            #pragma omp critical
            {
            freqs[node_idx] = counts_per_genomes;
            //}
            switch (std::get<0>(classified_kmers)) {
                case KmerType::CORE: {
                    std::for_each(
                        loop_collector.core_loops.begin(),
                        loop_collector.core_loops.end(),
                        [&](auto& core_kmer_freq_combs)
                        {
                        core_kmer_freq_combs[this->genome_number - 2] += 1;
                    });
                    break;
                }

                case KmerType::ACCESSORY: {
                    std::for_each(
                        loop_collector.accessory_loops.begin(),
                        loop_collector.accessory_loops.end(),
                        [&](auto& accessory_kmer_freq_combs)
                        {
                        accessory_kmer_freq_combs[this->genome_number - 2] += 1;
                     });
                    break;
                }

                case KmerType::UNIQUE: {
                    std::for_each(
                        loop_collector.unique_loops.begin(),
                        loop_collector.unique_loops.end(),
                        [&](auto& unique_kmer_freq_combs)
                        {
                        unique_kmer_freq_combs[this->genome_number - 2] += 1;
                     });
                    break;
                }
                case KmerType::NOT_PRESENT: {
                    break;
                }
            }
            }

        }
    }

    //     auto counts_per_genomes{this->count_kmer_per_genome(node_idx)};
    //     auto classified_kmers = this->classify_kmers(counts_per_genomes);
    //
    //     int reduced{std::reduce(counts_per_genomes.begin(), counts_per_genomes.end())};
    //
    //     if (reduced > 0) {
    //         freqs.push_back(counts_per_genomes);
    //     }
    //     //std::cout << "case:  " <<  kmer_type << std::endl;
    //     switch (std::get<0>(classified_kmers)) {
    //         case KmerType::CORE: {
    //             std::for_each(
    //                 loop_collector.core_loops.begin(),
    //                 loop_collector.core_loops.end(),
    //                 [&](auto& core_kmer_freq_combs)
    //                 {
    //                 core_kmer_freq_combs[this->genome_number - 2] += 1;
    //             });
    //             break;
    //         }
    //
    //         case KmerType::ACCESSORY: {
    //             std::for_each(
    //                 loop_collector.accessory_loops.begin(),
    //                 loop_collector.accessory_loops.end(),
    //                 [&](auto& accessory_kmer_freq_combs)
    //                 {
    //                 accessory_kmer_freq_combs[this->genome_number - 2] += 1;
    //              });
    //             break;
    //         }
    //
    //         case KmerType::UNIQUE: {
    //             std::for_each(
    //                 loop_collector.unique_loops.begin(),
    //                 loop_collector.unique_loops.end(),
    //                 [&](auto& unique_kmer_freq_combs)
    //                 {
    //                 unique_kmer_freq_combs[this->genome_number - 2] += 1;
    //              });
    //             break;
    //         }
    //         case KmerType::NOT_PRESENT: {
    //             break;
    //         }
    //     }
    // };
    return std::make_pair(freqs, std::ref(loop_collector));
};

