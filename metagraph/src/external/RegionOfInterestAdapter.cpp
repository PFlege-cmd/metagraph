//
// Created by Patrick Flege on 14/07/2026.
//

#include "RegionOfInterestAdapter.h"

#include "cli/graph_glue.hpp"
#include "graph/annotated_dbg.hpp"

struct HitsPerSequence {
    int* hitsPerSequence;
    int* sequenceOffsets;
    HitsPerSequence() = default;
    explicit HitsPerSequence(int* arr, int* offsets ) { hitsPerSequence = arr; sequenceOffsets = offsets; };
};
extern "C" {
    HitsPerSequence* retrieveKmersOfInterest(char * database) {
        auto results = new HitsPerSequence();
        graph_glue glue = graph_glue(0, NULL);
        static std::shared_ptr<mtg::graph::AnnotatedDBG> graph = glue.load_coord_dbg(database);

        return results;
    }


}