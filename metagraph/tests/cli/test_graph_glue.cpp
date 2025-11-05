//
// Created by Patrick Flege on 04/11/2025.
//

#include "test_graph_glue.h"

#include "gtest/gtest.h"
#include "mock_annotated_dbg.h"
#include "../../../../../../../opt/homebrew/Cellar/boost/1.88.0/include/boost/iostreams/filter/zstd.hpp"

TEST(pantools_interaction, call_graph) {
    mock_annotated_dbg mock_graph = mock_annotated_dbg();
    EXPECT_CALL(mock_graph, read_mapping_pantools_both_sides).Times(testing::AtMost(1));



}