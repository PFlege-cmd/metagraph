//
// Created by Patrick Flege on 20/11/2025.
//

#include "mock_graph_wrapper.h"


#include "cli/GraphWrapper.hpp"
#include "cli/SequenceRetriever.h"

#include "common/vector.hpp"
#include "tests/annotation/test_annotated_dbg_helpers.hpp"

#include <vector>
#include "gmock/gmock-function-mocker.h"
#include "gmock/gmock-spec-builders.h"


using testing::internal::MockFunction;
using testing::Return;
using testing::_;


TEST(TestGraphWrapper, testGetAnchorId) {
    mtg::graph::AnnotatedDBG* graph = nullptr;
    MockGraphWrapper wrapper(*graph);
    SequenceRetriever retriever(wrapper);

    const char * testsequence = "abcdefghijklmnopqrst";
    std::string anchor_sequence = testsequence;

    EXPECT_CALL(wrapper, retrieveAnchorId(testing::_)).WillOnce(Return(1L));
    retriever.retrieveAnchorId(anchor_sequence);

    ASSERT_EQ(true, true);
}

TEST(TestGraphWrapper, testGetAnchorCoords) {
    mtg::graph::AnnotatedDBG* graph = nullptr;
    MockGraphWrapper wrapper(*graph);
    SequenceRetriever retriever(wrapper);
    const char * testsequence = "abcdefghijklmnopqrst";
    std::string anchor_sequence = testsequence;
    const char * genome_name = "gene_the_gene";
    long long anchor_position = 2;

    retriever.setAnchorPosition(anchor_position);
    retriever.setGenome(genome_name);

    EXPECT_CALL(wrapper, retrieveAnchorCoordinates(testsequence, anchor_position, genome_name)).WillOnce(Return(anchor_position));

    retriever.retrieveAnchorId(anchor_sequence);
}

TEST(TestGraphWrapper, testGetFirstNodeOfSearchedSequenceId) {
    mtg::graph::AnnotatedDBG* graph = nullptr;
    MockGraphWrapper wrapper(*graph);
    SequenceRetriever retriever(wrapper);
    const char * testsequence = "abcdefghijklmnopqrst";
    std::string anchor_sequence = testsequence;
    const char * genome_name = "gene_the_gene";
    long long anchor_position = 2;
    long long sequenceStart = 10;
    retriever.setAnchorPosition(anchor_position);
    retriever.setGenome(genome_name);
    retriever.setSequenceStart(sequenceStart);

    uint64_t anchor_node_id = 21;
    uint64_t start_node_id = 42;

    EXPECT_CALL(wrapper, get_first_node_of_coord_range(anchor_node_id, anchor_position, std::string(genome_name), sequenceStart)).WillOnce(Return(start_node_id));

    uint64_t actual_start_id = retriever.retrieveStartId(anchor_node_id);

    ASSERT_EQ(start_node_id, actual_start_id);
}

TEST(TestGraphWrapper, testGetSequenceForRange) {
    mtg::graph::AnnotatedDBG* graph = nullptr;
    MockGraphWrapper wrapper(*graph);

    SequenceRetriever retriever(wrapper);
    const char * testsequence = "abcdefghijklmnopqrst";
    std::string anchor_sequence = testsequence;
    const char * genome_name = "gene_the_gene";
    long long anchor_position = 2;
    long long sequence_start= 10;
    long long sequence_end = 20;
    retriever.setAnchorPosition(anchor_position);
    retriever.setGenome(genome_name);
    retriever.setSequenceStart(sequence_start);
    retriever.setSequenceEnd(sequence_end);

    uint64_t start_node_id = 42;

    std::string result_region = "ACACACACAT";

    EXPECT_CALL(wrapper, get_sequence_for_coords(genome_name, sequence_start, sequence_end, start_node_id)).WillOnce(Return("ACACACACAT"));

    std::string retrieved_region = retriever.retrieveRegionForRange(start_node_id);

    ASSERT_EQ(result_region, retrieved_region);
}