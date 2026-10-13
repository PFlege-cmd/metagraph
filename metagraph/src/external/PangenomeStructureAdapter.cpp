//
// Created by Patrick Flege on 16/09/2026.
//

#include "PangenomeStructureAdapter.h"

#include "KmerClassifier.hpp"
#include "PangenomeStructureClassifier.h"
#include "cli/DeBruijnGraphWrapper.h"
#include "cli/GraphWrapper.hpp"
#include "cli/graph_glue.hpp"

#include <string>
namespace NativeStructures {
struct PangenomeStructureKmerMatrix;
}
extern "C" {
typedef PangenomeStructureAdapter::PangenomeStructureCollector kmerCollector;
/**
 * @brief Calculates kmer-based pangenome structure.
 *
 * This class runs, in parallel, the estimation of pangenome structure with increasing
 * numbers of genomes. It returns a specialized matrix having, for each loop iteration,
 * the number of core-, accessory-, or unique kmers. Each estimation into these categories
 * is done, for n genomes, in random combinations starting from n to 2 genomes. This
 * serves to later calculate pangenome-openness, based on Heap's Law
 *
 * @param database_path path to genome-database
 * @param loops number of times to loop over different random combinations of genomes
 * @return PangenomeStructureKmerMatrix for core-, accessory-, and unique kmers.
 */
NativeStructures::PangenomeStructureKmerMatrix* calculate_pangenome_structure(char * database_path, int loops) {
    std::string data_path = std::string(database_path);
    graph_glue glue = graph_glue(0, NULL);

    static std::shared_ptr<AnnotatedDBG> graph = glue.load_dbg(data_path);
    DeBruijnGraphWrapper graph_wrapper{*graph};
    int num_genomes{graph_wrapper.get_num_genomes()};
    std::cout << "num_genomes = " << num_genomes << std::endl;
    std::vector genome_labels{graph_wrapper.get_genome_labels()};
    unsigned long long number_nodes_primary_graph{graph_wrapper.get_number_nodes()/2 + 1};
    std::cout << "number_nodes_primary_graph = " << number_nodes_primary_graph << std::endl;

    std::vector<unsigned long long> all_nodes(graph_wrapper.get_number_nodes()/2);
    std::vector<long long> accessory_or_unique_node_ids{};
    std::iota(all_nodes.begin(), all_nodes.end(), 1);

    int core_threshold{num_genomes};
    int unique_threshold{1};
    KmerClassifier kmer_classifier{graph_wrapper, core_threshold, unique_threshold, num_genomes};

    kmerCollector collector;
    collector.core_loops = std::vector<std::vector<int>>(loops, std::vector<int>(num_genomes - 1, 0));
    collector.accessory_loops = std::vector<std::vector<int>>(loops, std::vector<int>(num_genomes - 1, 0));
    collector.unique_loops = std::vector<std::vector<int>>(loops, std::vector<int>(num_genomes - 1, 0));
    auto bit_freqs {boost::dynamic_bitset<>(num_genomes*all_nodes.size())};
    collector.freqs = std::move(bit_freqs);
    collector.number_nodes = all_nodes.size();
    //TODO: Modify Kmer-classfifier to avoid useless creation of vectors
    auto frequencies_and_collector{kmer_classifier.initialize_kmer_frequencies_bitvector(collector, all_nodes)};

    PangenomeStructureClassifier pangenome_classifier{graph_wrapper, num_genomes, loops};

    std::vector<std::vector<int>> combinations{pangenome_classifier.generate_genome_combinations_for_loop()};
    collector = pangenome_classifier.loop_for_increasing_pangenome_size(collector, combinations);

    auto kmerMatrix{PangenomeStructureAdapter::export_matrix(collector)};
    return kmerMatrix;
    }
}