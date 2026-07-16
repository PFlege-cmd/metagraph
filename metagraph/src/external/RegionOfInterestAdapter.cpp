//
// Created by Patrick Flege on 14/07/2026.
//

#include "RegionOfInterestAdapter.h"

#include "CoordinateRetriever.h"
#include "cli/graph_glue.hpp"
#include "graph/annotated_dbg.hpp"

struct HitsPerSequence {
    int* hitsPerSequence;
    int* sequenceOffsets;
    HitsPerSequence() = default;
    explicit HitsPerSequence(int* arr, int* offsets ) { hitsPerSequence = arr; sequenceOffsets = offsets; };
};
extern "C" {
    void retrieveKmersOfInterest(char ** kmers) {
        //char * database, int max_frequency, int number_kmers, char* kmers
        //auto results = new HitsPerSequence();
        //graph_glue glue = graph_glue(0, NULL);
        //static std::shared_ptr<mtg::graph::AnnotatedDBG> graph = glue.load_coord_dbg("database");

        //CoordinateRetriever coordinate_retriever(graph);
        std::cout << "Number of kmers:" << strlen(*kmers) << std::endl;
        //return results;
    }


}