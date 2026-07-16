//
// Created by Patrick Flege on 14/07/2026.
//

#ifndef GENOMECOORDINATE_H
#define GENOMECOORDINATE_H
#include "MatrixEntry.h"
#include "cli/GraphWrapper.hpp"


#include <string>


class GenomeCoordinate : public MatrixEntry {
public:
    typedef std::string Label;
    typedef std::vector<SmallVector<uint64_t>> genomeCoordinateVector;
    typedef std::tuple<Label, size_t, genomeCoordinateVector> genomeCoordinateTriple;

    GenomeCoordinate(std::string_view& genome, std::string_view& kmer, std::vector<int>& sequence_lengths, GraphWrapper& graph);
    ~GenomeCoordinate()= default;
    std::string_view& get_reverse_kmer() override;
    std::vector<std::array<int, 2>> locate() override;
    std::vector<std::array<int, 2>> locate_reverse() override;
    std::vector<std::array<int, 2>> locate_both_kmers() override;

private:
        std::string_view& reverse_kmer;
        std::string_view& create_kmer_reverse_complement(const std::string_view& kmer);
        void calculate_sequence_location(unsigned long coord, const std::vector<int>& sequence_lengths, std::array<int, 2>& position_and_location);
        void calculate_sequence_location_reverse(unsigned long coord, const std::vector<int>& sequence_lengths, std::array<int, 2>& position_and_location);

};



#endif //GENOMECOORDINATE_H
