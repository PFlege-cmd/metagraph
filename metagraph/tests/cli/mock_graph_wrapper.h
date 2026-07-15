//
// Created by Patrick Flege on 16/03/2026.
//

#pragma once
#include "cli/GraphWrapper.hpp"


#include "gmock/gmock-function-mocker.h"

using node_index = GraphWrapper::node_index_kmer;

class MockGraphWrapper : public GraphWrapper {
public:
    explicit MockGraphWrapper(mtg::graph::AnnotatedDBG& graph);
    explicit MockGraphWrapper(MockGraphWrapper& a);
  MOCK_METHOD(long,
              retrieveAnchorCoordinates,
              (std::string anchor_sequence, long anchor_position, std::string genome_name),
              (override));
    MOCK_METHOD(uint64_t, retrieveAnchorId, (std::string), (override));
    MOCK_METHOD((std::string), get_sequence_for_coords, (std::string, unsigned long long, unsigned long long, uint64_t start_index), (override));
    MOCK_METHOD(uint64_t, get_first_node_of_coord_range, (uint64_t, long long start_anchor, std::string, long long start_sequence), (override));
    MOCK_METHOD(unsigned long long, get_number_nodes, (), (override));
    MOCK_METHOD(kmer_frequencies, get_kmer_frequencies, (std::vector<node_index_kmer>), (override));
    MOCK_METHOD(int, get_num_genomes, (), (override, const));
    MOCK_METHOD(std::vector<GraphWrapper::genomeCoordinateTriple>, get_kmer_coordinates, (std::string), (override));
    MOCK_METHOD(size_t, get_K, (), (override, const));
};
