//
// Created by Patrick Flege on 14/07/2026.
//

#ifndef COORDINATERETRIEVER_H
#define COORDINATERETRIEVER_H
#include "MatrixEntry.h"
#include "MatrixEntryFactory.h"
#include "cli/GraphWrapper.hpp"


#include <string>


class CoordinateRetriever {
public:
    CoordinateRetriever(int no_of_sequences, int max_frequency, std::string& genome, int* sequences, GraphWrapper& graph);

    ~CoordinateRetriever();

    void set_graph(GraphWrapper& wrapper);
    int get_number_of_sequences();
    int get_sequence_length(int sequence_number);
    GraphWrapper& get_graph();
    int get_max_frequency();
    std::vector<std::array<int, 2>> get_kmer_positions(std::vector<std::unique_ptr<MatrixEntry>>&);
    int* get_sequences();
    //void setCoordinates(std::vector<std::unique_ptr<MatrixEntry>>& genomeCoordinates);
    std::vector<std::unique_ptr<MatrixEntry>>
    createEntries(MatrixEntryFactory& factory, const char** kmers, int number_of_kmers);
    std::vector<std::string_view> filter_max_frequency(std::vector<std::string_view> kmers);

private:
    int number_of_sequences;
    int max_frequency;
    std::string genome;
    int* sequences{};
    GraphWrapper& graph;
};



#endif //COORDINATERETRIEVER_H
