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
