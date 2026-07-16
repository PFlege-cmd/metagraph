//
// Created by Patrick Flege on 14/07/2026.
//

#include "CoordinateRetriever.h"
CoordinateRetriever::CoordinateRetriever(int number_of_sequences, int max_frequency, std::string& genome, int* sequences, std::vector<std::unique_ptr<MatrixEntry>>& genomeCoordinates, GraphWrapper& graph): number_of_sequences(number_of_sequences), max_frequency(max_frequency), genome(genome), sequences(sequences), graph(graph), genomeCoordinates(genomeCoordinates) {
    std::cout << "Created CoordinateRetriever" << std::endl;
    std::cout << "Max frequency: " << max_frequency << std::endl;
    std::cout << "Genome: " << genome << std::endl;
    std::cout << "Number of sequences: " << number_of_sequences << std::endl;
    std::cout << "Number of sequences in the graph: " << number_of_sequences << std::endl;
    std::cout << "Graph ref :"  << std::endl;
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
std::vector<std::array<int, 2>> CoordinateRetriever::get_kmer_positions() {
    std::vector<std::array<int, 2>> kmer_positions;
    std::for_each(this->genomeCoordinates.begin(), this->genomeCoordinates.end(), [&kmer_positions](auto const &genome_coords) {
        auto found_coords = genome_coords->locate_both_kmers();
        kmer_positions.insert(end(kmer_positions), begin(found_coords), end(found_coords));
    });
    return kmer_positions;
};
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

std::vector<std::unique_ptr<MatrixEntry>>& CoordinateRetriever::getCoordinates() {
    return genomeCoordinates;
}

std::vector<std::unique_ptr<MatrixEntry>>& CoordinateRetriever::createEmpty() {
    std::vector<std::unique_ptr<MatrixEntry>>* tmp
            = new std::vector<std::unique_ptr<MatrixEntry>>(0);
    return *tmp;
}