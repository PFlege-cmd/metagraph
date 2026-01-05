//
// Created by Patrick Flege on 04/11/2025.
//

#include "test_graph_glue.h"
#include "cli/graph_glue.hpp"
#include "cli/config/config.hpp"
#include "stdlib.h"
#include "gtest/gtest.h"
#include "mock_annotated_dbg.h"
#include "cli/cli_caller_interface.hpp"


using testing::MockFunction;
using testing::Return;
using testing::_;

//MockFunction<unique_ptr<Config>*> mock_function;
class MockCLICaller : public AbstractCommandLineInterface {
    public:
        MOCK_METHOD(unique_ptr<Config>, get_config, (int argc, char** cmd_arguments), ());
        MOCK_METHOD(int, call_flow, (const unique_ptr<Config>&), ());
};


TEST(pantools_interaction, get_cmd_args) {

    int argc = 5;
    char** argv = (char**)calloc(argc, sizeof(const char*));
    argv[0] = (char *) "/Users/patrick_flege/git/metagraph/metagraph/cmake-build-debug/metagraph_DNA5";
    argv[1] = (char*)"query";
    argv[2] = (char*)"--query-mode";
    argv[3] = (char*)"coords";
    argv[4] = (char*)"-i";
    graph_glue glue = graph_glue(argc, argv);
    MockCLICaller caller = MockCLICaller();
    AbstractCommandLineInterface* cli_caller = &caller;
    glue.set_cli_caller(*cli_caller);


    ASSERT_EQ(glue.get_cmd_arg_count(), 5);

}

TEST(pantool_interaction, create_config) {
    int argc = 4;
    char** argv = (char**)calloc(argc,  sizeof(const char*));
    argv[0] = (char *) "/Users/patrick_flege/git/metagraph/metagraph/cmake-build-debug/metagraph_DNA5";
    argv[1] = (char *) "query";
    argv[2] = (char *) "--query-mode";
    argv[3] = (char *) "coords";

    graph_glue glue = graph_glue(argc, argv);
    MockCLICaller caller = MockCLICaller();
    AbstractCommandLineInterface* cli_caller = &caller;
    glue.set_cli_caller(*cli_caller);

    EXPECT_CALL(caller, get_config(_, _)).Times(1);

    caller.get_config(argc, argv);

    ASSERT_EQ(1, 1);
}

TEST(pantool_interaction, call_cmdline_flow){
    int argc = 4;
    char** argv = (char**)calloc(argc,  sizeof(const char*));
    argv[0] = (char *) "/Users/patrick_flege/git/metagraph/metagraph/cmake-build-debug/metagraph_DNA5";
    argv[1] = (char *) "query";
    argv[2] = (char *) "--query-mode";
    argv[3] = (char *) "coords";

    graph_glue glue = graph_glue(argc, argv);
    MockCLICaller caller = MockCLICaller();
    AbstractCommandLineInterface* cli_caller = &caller;
    glue.set_cli_caller(*cli_caller);

    EXPECT_CALL(caller, call_flow(_)).Times(1);

    glue.call_cmdline_flow(nullptr);
}