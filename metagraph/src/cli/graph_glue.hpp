//
// Created by Patrick Flege on 23/07/2025.
//

#ifndef GRAPH_GLUE_H
#define GRAPH_GLUE_H
#include "graph/representation/succinct/dbg_succinct.hpp"


struct HitsPerSequence;
shared_ptr<mtg::graph::AnnotatedDBG> load_dbg();
void do_pantools_work(char* genome_name,
             int* sequence_lengths,
             int no_of_sequences,
             HitsPerSequence*& results,  char * read, std::shared_ptr<mtg::graph::AnnotatedDBG> graph
);
#endif //GRAPH_GLUE_H
