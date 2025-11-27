//
// Created by Patrick Flege on 20/11/2025.
//

//#ifndef GRAPHWRAPPER_H
//#define GRAPHWRAPPER_H
#pragma once

#include "common/vector.hpp"
#include "graph/annotated_dbg.hpp"

#include <vector>


//std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>>

class GraphWrapper {
    typedef uint64_t node_index;
    public:
        explicit GraphWrapper(mtg::graph::AnnotatedDBG* graph) {
            graph_ = graph;
        };
        virtual long retrieveAnchorCoordinates(std::string anchor_sequence,
                                                                     long anchor_position,
                                                                     std::string genome_name) = 0;
        virtual node_index retrieveAnchorId(std::string anchor_sequence) = 0;
        virtual const char * get_sequence_for_coords(std::string genome, unsigned long long start, unsigned long long end) = 0;
        virtual node_index get_first_node_of_coord_range(node_index anchor_index, long long start_anchor, std::string genome_name, long long start_sequence) = 0;
        mtg::graph::AnnotatedDBG * get_graph() const {
            return graph_;
        }

    protected:
        mtg::graph::AnnotatedDBG * graph_;
};



//#endif //GRAPHWRAPPER_H
