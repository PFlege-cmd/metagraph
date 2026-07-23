//
// Created by Patrick Flege on 14/07/2026.
//

#include "CoordinateRetriever.h"

#include "MatrixEntryFactory.h"
CoordinateRetriever::CoordinateRetriever(int number_of_sequences, int max_frequency, std::string& genome, int* sequences,GraphWrapper& graph): number_of_sequences(number_of_sequences), max_frequency(max_frequency), genome(genome), sequences(sequences), graph(graph) {

};

CoordinateRetriever::~CoordinateRetriever(){};
int CoordinateRetriever::get_number_of_sequences(){return this->number_of_sequences;};
int CoordinateRetriever::get_sequence_length(int sequence_number) {
    if (sequence_number < number_of_sequences && sequence_number >= 0) {
        return sequences[sequence_number];
    } else {
        throw std::invalid_argument("Invalid sequence_number");
    }
    return sequences[sequence_number];
};
std::vector<std::array<int, 2>> CoordinateRetriever::get_kmer_positions(std::vector<std::unique_ptr<MatrixEntry>>& genomeCoordinates) {
    std::vector<std::array<int, 2>> kmer_positions;
    std::for_each(genomeCoordinates.begin(), genomeCoordinates.end(), [&kmer_positions](auto const &genome_coords) {
        auto found_coords = genome_coords->locate_both_kmers();
        kmer_positions.insert(end(kmer_positions), begin(found_coords), end(found_coords));
    });
    return kmer_positions;
};

void CoordinateRetriever::set_graph(GraphWrapper& wrapper) {
    this->graph = wrapper;
}

GraphWrapper& CoordinateRetriever::get_graph() {
    return graph;
}
int CoordinateRetriever::get_max_frequency() {
    return max_frequency;
}

int* CoordinateRetriever::get_sequences() {
    return sequences;
}

// void CoordinateRetriever::setCoordinates(
//         std::vector<std::unique_ptr<MatrixEntry>>& genomeCoordinates) {
//     // std::for_each(genomeCoordinates.begin(), genomeCoordinates.end(), [&](auto &genome_coords) {
//     //     this->genomeCoordinates.push_back(std::move(*genome_coords));
//     // });
//
//     this->genomeCoordinates.insert(this->genomeCoordinates.end(), genomeCoordinates.begin(), genomeCoordinates.end());
// }

std::vector<std::unique_ptr<MatrixEntry>>
CoordinateRetriever::createEntries(MatrixEntryFactory& factpry,
                                   const char** kmers,
                                   int number_of_kmers) {
    return factpry.createKmerEntries(kmers, number_of_kmers);
}

std::vector<std::string_view> CoordinateRetriever::filter_max_frequency(std::vector<std::string_view> kmers) {
    std::vector<std::string_view> filtered = std::vector<std::string_view>();
    std::for_each(kmers.begin(), kmers.end(), [&](auto kmer) {
        auto freqs = this->graph.get_kmer_frequencies(kmer); // OK since we only get a single entry

        for (size_t i = 0; i < freqs.size(); i++) {
            auto lab = std::get<0>(freqs[i]);
            auto equal_to_genome = this->genome.compare(lab);
            if (equal_to_genome == 0) {
                if (static_cast<int>(std::get<2>(freqs[i])[0]) <= this->max_frequency) {
                    filtered.push_back(kmer);
                }
            }
        }
    });
    return filtered;
}