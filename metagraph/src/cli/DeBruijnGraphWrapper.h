//
// Created by Patrick Flege on 26/11/2025.
//

#ifndef DEBRUIJNGRAPHWRAPPER_H
#define DEBRUIJNGRAPHWRAPPER_H
#include "GraphWrapper.hpp"


class DeBruijnGraphWrapper : public GraphWrapper {
    public:
    typedef std::string Label;
    typedef uint64_t node_index;
    explicit DeBruijnGraphWrapper(mtg::graph::AnnotatedDBG &graph);
    ~DeBruijnGraphWrapper();
    long retrieveAnchorCoordinates(std::string anchor_sequence,
                                   long anchor_position,
                                   std::string genome_name) override;
    uint64_t retrieveAnchorId(std::string anchor_sequence) override;
    std::string get_sequence_for_coords(std::string genome,
                                        unsigned long long start,
                                        unsigned long long end,
                                        node_index start_index) override;
    uint64_t get_first_node_of_coord_range(uint64_t anchor_index, long long start_anchor, std::string genome_name, long long start_sequence) override;

    unsigned long long get_number_nodes() override;

    kmer_frequencies get_kmer_frequencies(std::vector<node_index> nodes) override;
    int get_num_genomes() const override;
};



#endif //DEBRUIJNGRAPHWRAPPER_H
