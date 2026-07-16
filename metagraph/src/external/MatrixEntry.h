//
// Created by Patrick Flege on 15/07/2026.
//

#ifndef MATRIXENTRY_H
#define MATRIXENTRY_H


#include "MatrixEntry.h"

#include "cli/GraphWrapper.hpp"

#include <vector>
class MatrixEntry {
public:
    typedef std::string Label;
    typedef std::vector<SmallVector<uint64_t>> genomeCoordinateVector;
    typedef std::tuple<Label, size_t, genomeCoordinateVector> genomeCoordinateTriple;

    explicit MatrixEntry(std::string_view& genome, std::string_view& kmer, std::vector<int>& sequence_lengths, GraphWrapper& graph)
    :genome(genome),
    kmer(kmer),
    sequence_lengths(sequence_lengths),
    graph(graph) {};

    virtual ~MatrixEntry() = default;
    virtual std::string_view& get_reverse_kmer() = 0;
    virtual std::vector<std::array<int, 2>> locate() = 0;
    virtual std::vector<std::array<int, 2>> locate_reverse() = 0;
    virtual std::vector<std::array<int, 2>> locate_both_kmers() = 0;

protected:
    std::string_view& genome;
    std::string_view& kmer;
    std::vector<int>& sequence_lengths;
    GraphWrapper& graph;
};


#endif //MATRIXENTRY_H
