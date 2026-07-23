#include "mock_matrix_entry.h"
//
// Created by Patrick Flege on 16/07/2026.
//
 MockMatrixEntry::MockMatrixEntry(MockMatrixEntry& other) : MatrixEntry(other.genome, other.kmer, other.sequence_lengths, other.graph) {
    /// explicit MockMatrixEntry(std::string_view& genome, std::string_view& kmer,
    /// std::vector<int>& sequence_lengths, GraphWrapper& graph);
    //explicit MockMatrixEntry(MockMatrixEntry & other);
     this->graph = other.graph;
     this->genome = other.genome;
     this->kmer = other.kmer;
     this->sequence_lengths = other.sequence_lengths;
 }

 MockMatrixEntry::MockMatrixEntry(std::string_view& genome,
                                  std::string_view& kmer,
                                  std::vector<int>& sequence_lengths,
                                  GraphWrapper& graph) : MatrixEntry(genome, kmer,sequence_lengths, graph){
 }

