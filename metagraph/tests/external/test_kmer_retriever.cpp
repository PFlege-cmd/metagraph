#include "asio/error.hpp"
#include "external/CoordinateRetriever.h"


#include "gtest/gtest.h"
#include "tests/cli/mock_graph_wrapper.h"

//
// Created by Patrick Flege on 14/07/2026.
//
struct HitsPerSequence {
    int* hitsPerSequence;
    int* sequenceOffsets;
    HitsPerSequence() = default;
    explicit HitsPerSequence(int* arr, int* offsets ) { hitsPerSequence = arr; sequenceOffsets = offsets; };
};
TEST(testRetrieveKmers, TestRegionOfInterest) {
    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    int max_frequency = 10;
    std::string genome = "genome1";
    int* seqs = new int[3] {10, 20, 50 };
    CoordinateRetriever* retriever
            = new CoordinateRetriever(3, max_frequency, genome, seqs, wrapper);
    ASSERT_EQ(3, retriever->get_number_of_sequences());
}

TEST(testGetSequenceLength, TestRegionOfInterest) {
    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    int max_frequency = 10;
    std::string genome = "genome1";
    int* seqs = new int[3] {10, 20, 50 };
    CoordinateRetriever* retriever
            = new CoordinateRetriever(3, max_frequency, genome, seqs, wrapper);

    ASSERT_EQ(10, retriever->get_sequence_length(0));
    ASSERT_EQ(20, retriever->get_sequence_length(1));
    ASSERT_EQ(50, retriever->get_sequence_length(2));
    ASSERT_THROW(retriever->get_sequence_length(-1), std::invalid_argument);
    ASSERT_THROW(retriever->get_sequence_length(3), std::invalid_argument);
}



