//
// Created by Patrick Flege on 23/07/2025.
//

#ifndef GRAPH_GLUE_H
#define GRAPH_GLUE_H
#include "cli_caller.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/annotated_dbg.hpp"
#include "cli/config/config.hpp"

using namespace mtg::graph;
using namespace mtg::cli;
struct HitsPerSequence;
class graph_glue {
    public:
    graph_glue(int argc, char** argv);
    ~graph_glue();
    //struct HitsPerSequence;
    shared_ptr<AnnotatedDBG> load_dbg();
    void do_pantools_work(char* genome_name,
                 int* sequence_lengths,
                 int no_of_sequences,
                 HitsPerSequence*& results,  const char * read, const std::shared_ptr<AnnotatedDBG>& graph
    );
    char** get_cmd_arguments();
    int get_cmd_arg_count();
    unique_ptr<Config> create_config();
    void call_cmdline_flow(const std::unique_ptr<mtg::cli::Config>& config);
    void set_cmd_arguments(char* arguments[]);
    void set_cli_caller(AbstractCommandLineInterface& cli);

    private:
        char** cmd_arguments;
        int argc;
        AbstractCommandLineInterface* cli_caller;
};

struct HitsPerSequence;
//struct HitsPerSequence;
//shared_ptr<mtg::graph::AnnotatedDBG> load_dbg();
//void do_pantools_work(char* genome_name,
//             int* sequence_lengths,
//             int no_of_sequences,
//             HitsPerSequence*& results,  char * read, std::shared_ptr<mtg::graph::AnnotatedDBG> graph
//);
#endif //GRAPH_GLUE_H
