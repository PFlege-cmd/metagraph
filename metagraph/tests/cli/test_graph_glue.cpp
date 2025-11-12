//
// Created by Patrick Flege on 04/11/2025.
//

#include "test_graph_glue.h"
#include "cli/graph_glue.hpp"

#include "gtest/gtest.h"
#include "mock_annotated_dbg.h"
#include "../../../../../../../opt/homebrew/Cellar/boost/1.88.0/include/boost/iostreams/filter/zstd.hpp"

TEST(pantools_interaction, get_cmd_args) {

    int argc = 5;
    char** argv = (char**)malloc(argc * sizeof(const char*));
    argv[0] = (char *) "/Users/patrick_flege/git/metagraph/metagraph/cmake-build-debug/metagraph_DNA5";
    argv[1] = (char*)"query";
    argv[2] = (char*)"--query-mode";
    argv[3] = (char*)"coords";
    argv[4] = (char*)"-i";
    graph_glue glue = graph_glue(argc, argv);


    ASSERT_EQ(glue.get_cmd_arg_count(), 5);

}