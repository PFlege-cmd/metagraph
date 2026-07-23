//
// Created by Patrick Flege on 20/11/2025.
//

#pragma once

#include "common/vector.hpp"
#include "graph/annotated_dbg.hpp"

#include <vector>

class GraphWrapper {
    public:
    typedef std::string Label;
    typedef std::vector<SmallVector<uint64_t>> genomeCoordinateVector;
    typedef std::tuple<Label, size_t, genomeCoordinateVector> genomeCoordinateTriple;
    typedef uint64_t node_index_kmer;
    typedef std::vector<std::tuple<mtg::graph::AnnotatedSequenceGraph::Label, size_t, std::vector<size_t>>> kmer_frequencies;

    virtual ~GraphWrapper() = default;
    explicit GraphWrapper(mtg::graph::AnnotatedDBG* graph) {
            graph_ = graph;
        };
        virtual long retrieveAnchorCoordinates(std::string anchor_sequence,
                                                                     long anchor_position,
                                                                     std::string genome_name) = 0;
        virtual node_index_kmer retrieveAnchorId(std::string anchor_sequence) = 0;
        virtual std::string get_sequence_for_coords(std::string genome,
                                                    unsigned long long start,
                                                    unsigned long long end,
                                                    node_index_kmer start_index) = 0;
        virtual node_index_kmer get_first_node_of_coord_range(node_index_kmer anchor_index, long long start_anchor, std::string genome_name, long long start_sequence) = 0;
        mtg::graph::AnnotatedDBG * get_graph() const {
            return graph_;
        }
        virtual unsigned long long get_number_nodes() = 0;
        virtual kmer_frequencies get_kmer_frequencies(std::vector<node_index_kmer> nodes) = 0;
        virtual kmer_frequencies get_kmer_frequencies(std::string_view kmers) = 0;

        virtual int get_num_genomes() const = 0;
        virtual std::vector<genomeCoordinateTriple> get_kmer_coordinates(std::string_view& kmer) = 0;
        virtual size_t get_K() const = 0;
    protected:
        mtg::graph::AnnotatedDBG * graph_;
};