//
// Created by Patrick Flege on 14/07/2026.
//

#include "GenomeCoordinate.h"

#include "../cli/GraphWrapper.hpp"
#include "common/seq_tools/reverse_complement.hpp"

GenomeCoordinate::GenomeCoordinate(std::string& genome,
                                   std::string& kmer,
                                   std::vector<int>& sequence_lengths,
                                   GraphWrapper& graph):genome(genome), kmer(kmer), reverse_kmer(create_kmer_reverse_complement(kmer)), sequence_lengths(sequence_lengths), graph(graph) {
}

std::string& GenomeCoordinate::get_reverse_kmer() {
    return reverse_kmer;
}

std::string& GenomeCoordinate::create_kmer_reverse_complement(const std::string& kmer) {
    char reverse_kmer[strlen(kmer.data())];
    memcpy(reverse_kmer, kmer.data(), kmer.size());
    auto rev_kmer = new std::string(reverse_kmer);
    reverse_complement(rev_kmer->begin(), rev_kmer->end());
    return *rev_kmer;
};


std::vector<std::array<int, 2>> GenomeCoordinate::locate() {
    std::string current_kmer = kmer;
    std::vector<std::array<int, 2>> kmer_positions;

        std::vector<genomeCoordinateTriple> coordinates = this->graph.get_kmer_coordinates(current_kmer);
        for (unsigned long coord_idx = 0; coord_idx < coordinates.size();coord_idx++) {
            std::string current_genome = std::get<0>(coordinates[coord_idx]);
            if (current_genome != genome) {
                continue;
            }
            auto coords_for_genome = std::get<2>(coordinates[coord_idx]);
            auto test = coords_for_genome[0];
            auto coords_for_genome_size = coords_for_genome[0].size();
            for (unsigned long current_coord_idx = 0; current_coord_idx < coords_for_genome_size; current_coord_idx++) {
                std::array<int,2>  target_and_position = {0, 0};
                calculate_sequence_location(coords_for_genome[0][current_coord_idx],
                                            sequence_lengths, target_and_position);
                int target_sequence = target_and_position[1]; // TODO: Make it inline with pantools
                int pantools_location = target_and_position[0];
                int loc = pantools_location;

                if (loc >= 0 && loc <= sequence_lengths[target_sequence]) {
                    kmer_positions.push_back({target_sequence + 1, loc});
                }
            }
        }
    return kmer_positions;
}

void GenomeCoordinate::calculate_sequence_location(unsigned long coord, const std::vector<int>& sequence_lengths, std::array<int, 2>& position_and_location) {
    int current_length = 0;
    int k = this->graph.get_K();
    int target_sequence = -1;
    for (int i = 0; i < (int)sequence_lengths.size(); i++) {
        int seq_length = sequence_lengths[i] - k + 1;
        if ((int) coord >= current_length && (int) coord < current_length + seq_length) {
            target_sequence = i;
            break;
        }
        current_length += seq_length;
    }

    int pantools_location = coord - current_length;
    position_and_location[0] = pantools_location;
    position_and_location[1] = target_sequence;
}

void GenomeCoordinate::calculate_sequence_location_reverse(unsigned long coord, const std::vector<int>& sequence_lengths, std::array<int, 2>& position_and_location) {
    int current_length = 0;
    int k = this->graph.get_graph()->get_graph().get_k();
    int target_sequence = -1;
    for (int i = 0; i < (int)sequence_lengths.size(); i++) {
        int seq_length = sequence_lengths[i] - k + 1;
        if ((int) coord >= current_length && (int) coord < current_length + seq_length) {
            target_sequence = i;
            break;
        }
        current_length += seq_length;
    }

    int pantools_location = coord - current_length;
    position_and_location[0] = pantools_location;
    position_and_location[1] = target_sequence;
}
