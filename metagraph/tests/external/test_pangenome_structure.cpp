#include "cli/DeBruijnGraphWrapper.h"
#include "cli/graph_glue.hpp"
#include "external/KmerClassifier.hpp"
#include "external/KmerUtils.h"
#include "external/PangenomeStructureAdapter.h"
#include "external/PangenomeStructureClassifier.h"
#include "tests/cli/mock_graph_wrapper.h"

#include "gtest/gtest.h"
//
// Created by Patrick Flege on 16/09/2026.
//

typedef PangenomeStructureAdapter::PangenomeStructureCollector kmerCollector;
TEST(testPangenomeStructure, testIntegration) {

    int loops{100};
    std::string data_path = std::string("../../../patrick-pan-tools/covid_DB");
    graph_glue glue = graph_glue(0, NULL);

    static std::shared_ptr<AnnotatedDBG> graph = glue.load_dbg(data_path);
    DeBruijnGraphWrapper graph_wrapper{*graph};
    int num_genomes{graph_wrapper.get_num_genomes()};

    std::cout << "num_genomes = " << num_genomes << std::endl;
    std::vector genome_labels{graph_wrapper.get_genome_labels()};
    unsigned long long number_nodes_primary_graph{graph_wrapper.get_number_nodes()/2 + 2};
    std::cout << "number_nodes_primary_graph = " << number_nodes_primary_graph << std::endl;

    std::vector<unsigned long long> all_nodes(graph_wrapper.get_number_nodes()/2);
    std::iota(all_nodes.begin(), all_nodes.end(), 1);

    int core_threshold{num_genomes};
    int unique_threshold{1};
    KmerClassifier kmer_classifier{graph_wrapper, core_threshold, unique_threshold, num_genomes};

    kmerCollector collector;
    collector.core_loops = std::vector<std::vector<int>>(loops, std::vector<int>(num_genomes - 1, -1));
    collector.accessory_loops = std::vector<std::vector<int>>(loops, std::vector<int>(num_genomes - 1, -1));
    collector.unique_loops = std::vector<std::vector<int>>(loops, std::vector<int>(num_genomes - 1, -1));

    auto accessory_unique_kmers_and_collector{kmer_classifier.initialize_core_kmers(collector, all_nodes)};

    std::vector<std::vector<int>> accessory_unique_kmers{accessory_unique_kmers_and_collector.first};
    constexpr int unique_kmers = 64;
    constexpr int accessory_kmers = 247;
    constexpr int core_kmers = 89282;

    for (size_t i = 0; i < static_cast<size_t> (loops); i++) {
        ASSERT_EQ(accessory_unique_kmers_and_collector.second.core_loops[i][num_genomes - 2], core_kmers);
        ASSERT_EQ(accessory_unique_kmers_and_collector.second.accessory_loops[i][num_genomes - 2], accessory_kmers);
        ASSERT_EQ(accessory_unique_kmers_and_collector.second.unique_loops[i][num_genomes - 2], unique_kmers);
    }
}

TEST(testPangenomeStructure, testGenerateRandomCombinations) {
    int loops{100};
    int number_genomes{4};

    AnnotatedDBG* graph = nullptr;
    MockGraphWrapper wrapper = MockGraphWrapper(*graph);

    PangenomeStructureClassifier pangenome_classifier{wrapper, number_genomes, loops};

    std::mt19937 random_engine{std::random_device{}()};
    std::vector<std::uniform_int_distribution<int>> distributions(number_genomes);
    for (size_t i = 0; i < static_cast<size_t>(number_genomes); i++) {
        distributions[i] = std::uniform_int_distribution<int>{0, static_cast<int>(i)};
    }

    auto genome_combinations{pangenome_classifier.generate_genome_combinations(random_engine, distributions)};
    for (size_t i = 0; i < genome_combinations.size(); i++) {
        for (size_t j = 0; j < genome_combinations.size() && i != j; j++) {
            ASSERT_NE(genome_combinations[j], genome_combinations[i]);
        }
    }
}

TEST(testPangenomeStructure, testGenerateRandomCombinationsforLoops) {
    int loops{100};
    int number_genomes{4};

    AnnotatedDBG* graph = nullptr;
    MockGraphWrapper wrapper = MockGraphWrapper(*graph);

    PangenomeStructureClassifier pangenome_classifier{wrapper, number_genomes, loops};

    std::vector<std::vector<int>> combinations{pangenome_classifier.generate_genome_combinations_for_loop()};
    ASSERT_TRUE(!combinations.empty());
    ASSERT_EQ(combinations.size(), loops);
    for (size_t i = 0; i < combinations.size(); i++) {
        for (size_t j = 0; j < combinations[i].size(); j++) {
            for (size_t k = 0; k < combinations[i].size() && k != j; k++) {
                ASSERT_NE(combinations[i][k], combinations[i][j]);
            }
        }
    }
}

TEST(testPangenomeStructure, testClassifyKmersForLimitedGenomes) {

    std::vector<int> genome_idx{1, 2, 3};

    std::vector<int> kmer_frequencies_accessory{2, 0, 2, 3};
    std::vector<int> kmer_frequencies_core{0, 1, 2, 3};
    auto kmer_classified{KmerUtils::classify_selected_kmers(kmer_frequencies_accessory, genome_idx)};
}