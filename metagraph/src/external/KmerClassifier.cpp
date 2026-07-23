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
        //int highest =  outer_genome_idx.second > inner_genome_idx.second?outer_genome_idx.second : inner_genome_idx.second;

        if (lowest > 0) {
            ++distinct_shared[outer_index][inner_index];
        }
        //++this->distinct_total[outer_index][inner_index];

        this->all_shared_matrix[outer_index][inner_index] += lowest;
        //this->all_total_matrix[outer_index][inner_index] += highest;
    });
}

int* KmerClassifier::flatten(const std::vector<std::vector<int>>& matrix) {
    // std::vector<std::vector<const int>::iterator> all_shared_iter_start;
    //
    // all_shared_iter_start.resize(genome_number);
    // int* flattened = new int[genome_number*genome_number];
    //
    // for (int i = 0; i < genome_number; i++) {
    //     all_shared_iter_start[i] = matrix[i].begin();
    // }
    //
    // for (int i = 0; i < genome_number; i++) {
    //     for (int j = 0; j < genome_number; j++) {
    //         std::cout  << *all_shared_iter_start[i] << ","  << std::endl;
    //         flattened[i*genome_number + j] = *all_shared_iter_start[i];
    //         ++all_shared_iter_start[i];
    //     }
    // }

    int* flattened = new int[genome_number*genome_number];

    for (int i = 0; i < genome_number; i++) {
        for (int j = 0; j < genome_number; j++) {
            flattened[i*genome_number + j] = matrix[i][j];
        }
    }
    return flattened;
}


//Note: this needs to call the num_genomes beforehand.
// template <int genome_number>
// std::array<int, genome_number> KmerClassifier::classify_kmer(node_index_kmer idx) {
//     std::array<int, genome_number> kmer_counts(0);
//     auto frequencies = this->get_kmer_frequencies(std::vector<node_index_kmer>(idx));
//     for (size_t i = 0; i < frequencies.size(); i++) {
//         auto genome_name =  std::get<0>(frequencies[i]);
//         std::regex contig_del("contig");
//         std::sregex_token_iterator it(genome_name.begin(), genome_name.end(), contig_del, -1);
//         std::sregex_token_iterator end;
//
//         //int genome_number = std::get<0>(frequencies[i]).
//     }
//     return kmer_counts;
// }
