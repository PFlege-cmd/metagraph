//
// Created by Patrick Flege on 20/11/2025.
//

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

class MockGraphWrapper : public GraphWrapper {
    public:
        explicit MockGraphWrapper(mtg::graph::AnnotatedDBG& graph);
        MOCK_METHOD(long,retrieveAnchorCoordinates, (std::string anchor_sequence, long anchor_position, std::string genome_name), (override));
        MOCK_METHOD(uint64_t, retrieveAnchorId, (std::string), (override));
        MOCK_METHOD((const char *), get_sequence_for_coords, (std::string, unsigned long long, unsigned long long), (override));
        MOCK_METHOD(uint64_t, get_first_node_of_coord_range, (uint64_t, long long start), (override));
        //std::vector<std::tuple<std::string, size_t, std::vector<SmallVector<uint64_t>>>> retrieveAnchorCoordinates(char* anchor_sequence, long anchor_position, char* genome_name) override;
};

MockGraphWrapper::MockGraphWrapper(mtg::graph::AnnotatedDBG &graph) : GraphWrapper(&graph) {
    //I think that &stuff is converted, ie copied, by copying tyhr pointer to it?

    std::cout << "MockGraphWrapper::MockGraphWrapper()" << &graph  << std::endl;
}


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

    EXPECT_CALL(wrapper, get_first_node_of_coord_range(anchor_node_id, sequenceStart)).WillOnce(Return(start_node_id));

    uint64_t actual_start_id = retriever.retrieveStartId(anchor_node_id);

    ASSERT_EQ(start_node_id, actual_start_id);
}