//
// Created by Patrick Flege on 17/07/2026.
//

#include "GenomeCoordinateFactory.h"

#include "GenomeCoordinate.h"
std::unique_ptr<MatrixEntryFactory> GenomeCoordinateFactory::instance = nullptr;

GenomeCoordinateFactory::GenomeCoordinateFactory(std::string genome,
                                                 int no_seqs,
                                                 int* seqs,
                                                 GraphWrapper& wrapper)
    : wrapper(wrapper) {
    this->seqs = std::vector<int>(no_seqs); // Needs to be initialized, else BAD ACCESS
    this->genome = genome;
    for (int i = 0; i < no_seqs; i++) {

        this->seqs[i] = seqs[i];
    }
};

std::vector<std::unique_ptr<MatrixEntry>> GenomeCoordinateFactory::createKmerEntries(const char** kmers, int numberOfKmers) {
    std::cout << "Creating kmer entries..." << std::endl;
    std::vector<int> idxs = std::vector<int>(numberOfKmers, 0);
    std::vector<std::unique_ptr<MatrixEntry>> coordinates;
    std::iota(idxs.begin(),  idxs.end(), 0); //TODO:: GOes bad here...
    std::for_each(idxs.begin(),  idxs.end(), [&](int idx) {
        std::cout << idx <<std::endl;
        std::cout << kmers[idx] << std::endl;
            auto kmer_view = std::string_view(kmers[idx]);
            std::cout << "Kmerview is:" << std::endl;

            std::cout << kmer_view << std::endl;
            auto coordinate = std::make_unique<GenomeCoordinate>(this->genome, kmer_view, this->seqs, this->wrapper);
            std::cout << "Coordinate is:" << coordinate->get_kmer() <<  std::endl;
            coordinates.push_back(std::move(coordinate));
            //coordinates.push_back(std::make_unique<GenomeCoordinate>(this->genome, kmer_view, this->seqs, this->wrapper));

    });
    for (size_t i = 0; i < coordinates.size(); i++) {
        std::cout << coordinates[i].get()->get_reverse_kmer() << std::endl;
        std::cout << coordinates[i].get()->get_kmer() << std::endl;

    }
    std::cout << "Return" << std::endl;
    return coordinates;
}


std::unique_ptr<MatrixEntryFactory>& GenomeCoordinateFactory::getInstance(std::string genome,int no_seqs, int* seqs, GraphWrapper& wrapper) {
    if (!GenomeCoordinateFactory::instance) {
        GenomeCoordinateFactory::instance = std::make_unique<GenomeCoordinateFactory>(genome, no_seqs, seqs, wrapper);
    }
    return GenomeCoordinateFactory::instance;
}