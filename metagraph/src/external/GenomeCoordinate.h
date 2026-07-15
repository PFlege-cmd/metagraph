//
// Created by Patrick Flege on 14/07/2026.
//

#ifndef GENOMECOORDINATE_H
#define GENOMECOORDINATE_H
#include "cli/GraphWrapper.hpp"


#include <string>


class GenomeCoordinate {
public:
    typedef std::string Label;
    typedef std::vector<SmallVector<uint64_t>> genomeCoordinateVector;
    typedef std::tuple<Label, size_t, genomeCoordinateVector> genomeCoordinateTriple;

    GenomeCoordinate(std::string& genome, std::string& kmer, std::vector<int>& sequence_lengths, GraphWrapper& graph);
    ~GenomeCoordinate()= default;
    std::string& get_reverse_kmer();
    std::vector<std::array<int, 2>> locate();

private:
        std::string& genome;
        std::string& kmer;
        std::string& reverse_kmer;
        std::vector<int>& sequence_lengths;
        GraphWrapper& graph;
        std::string& create_kmer_reverse_complement(const std::string& kmer);
        void calculate_sequence_location(unsigned long coord, const std::vector<int>& sequence_lengths, std::array<int, 2>& position_and_location);
        void calculate_sequence_location_reverse(unsigned long coord, const std::vector<int>& sequence_lengths, std::array<int, 2>& position_and_location);

};



#endif //GENOMECOORDINATE_H
