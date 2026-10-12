#include "../../src/external/CombinationUtils.h"
#include "cli/DeBruijnGraphWrapper.h"
#include "cli/graph_glue.hpp"
#include "external/KmerClassifier.hpp"
#include "external/KmerUtils.h"
#include "external/PangenomeStructureAdapter.h"
#include "external/PangenomeStructureClassifier.h"
#include "tests/cli/mock_graph_wrapper.h"

#include "gtest/gtest.h"
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
//
// Created by Patrick Flege on 16/09/2026.
//

vector<vector<int>> create_five_genome_test_combinations() {
    std::vector<int> combinations1 { 0, 1, 2, 3 };
    std::vector<int> combinations2 { 0, 3, 2, 1 };
    std::vector<int> combinations3 { 2, 1, 0, 3 };
    std::vector<int> combinations4 { 3, 2, 1, 1 };
    std::vector<int> combinations5 { 0, 1, 2, 3 };
    vector<vector<int>> all_combinations { combinations1, combinations2, combinations3,
                                           combinations4, combinations5 };
    return all_combinations;
}

std::vector<std::vector<int>> create_test_frequencies_five_genomes(size_t frequencies) {
    std::mt19937 random_engine { std::random_device {}() };
    std::uniform_int_distribution<int> distribution { 0, 10 };

    std::vector<std::vector<int>> freqs(frequencies, std::vector<int>(5, 0));

    for (size_t i = 0; i < frequencies; i++) {

        freqs[i][0] = distribution(random_engine);
        freqs[i][1] = distribution(random_engine);
        freqs[i][2] = distribution(random_engine);
        freqs[i][3] = distribution(random_engine);
        freqs[i][4] = distribution(random_engine);

    }
    return freqs;
}

boost::dynamic_bitset<>& fill_bitset_w_freqs(
        std::vector<std::vector<int>> freqs,
        boost::dynamic_bitset<>& freqs_bitset) {


    for (size_t f_1 = 0; f_1 < freqs.size(); f_1++) {
        for (size_t f_2 = 0; f_2 < 5; f_2++) {
            if (freqs[f_1][f_2] != 0) {
                freqs_bitset.flip(f_1 * 5 + f_2);
            }
        }
    }
    return freqs_bitset;
}
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

    boost::dynamic_bitset<> freqs_bitset = boost::dynamic_bitset<>(all_nodes.size()*num_genomes);
    boost::dynamic_bitset<>& bitset_ref{freqs_bitset};

    collector.freqs = bitset_ref;
    collector.number_nodes = all_nodes.size();

    auto accessory_unique_kmers_and_collector{kmer_classifier.initialize_kmer_frequencies(collector, all_nodes)};

    constexpr int unique_kmers = 64;
    constexpr int accessory_kmers = 405;
    constexpr int core_kmers = 29700;

    for (size_t i = 0; i < static_cast<size_t> (loops); i++) {
        ASSERT_EQ(accessory_unique_kmers_and_collector.core_loops[i][num_genomes - 2], core_kmers);
        ASSERT_EQ(accessory_unique_kmers_and_collector.accessory_loops[i][num_genomes - 2], accessory_kmers);
        ASSERT_EQ(accessory_unique_kmers_and_collector.unique_loops[i][num_genomes - 2], unique_kmers);
    }
    PangenomeStructureClassifier pangenome_classifier{graph_wrapper, num_genomes, loops};

    std::vector<std::vector<int>> combinations{pangenome_classifier.generate_genome_combinations_for_loop()};

    ASSERT_EQ(combinations.size(), loops);

    collector = pangenome_classifier.loop_for_increasing_pangenome_size(collector, combinations);

    for (int loop = 0; loop < loops; loop++) {
        ASSERT_GT(std::reduce(collector.accessory_loops[loop].begin(), collector.accessory_loops[loop].end()), 0);
    }

    auto kmerMatrix{*PangenomeStructureAdapter::export_matrix(collector)};
    ASSERT_NE(kmerMatrix.core, nullptr);
    ASSERT_NE(kmerMatrix.accessory, nullptr);
    ASSERT_NE(kmerMatrix.unique, nullptr);
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

    auto genome_combinations{CombinationUtils::generate_combinations(random_engine, distributions, number_genomes)};
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

    std::vector<int> kmer_frequencies_accessory{2, 0, 2, 10};
    std::vector<int> kmer_frequencies_core{0, 14, 2, 11};
    auto accessory_kmer_classified{KmerUtils::classify_selected_kmers(kmer_frequencies_accessory, genome_idx)};
    auto core_kmer_classified{KmerUtils::classify_selected_kmers(kmer_frequencies_core, genome_idx)};

    ASSERT_EQ(std::get<0>(accessory_kmer_classified), KmerType::ACCESSORY);
    ASSERT_EQ(std::get<0>(core_kmer_classified), KmerType::CORE);

    auto accessory_pair_vector{std::get<1>(accessory_kmer_classified)};
    auto core_pair_vector{std::get<1>(core_kmer_classified)};

    ASSERT_EQ(accessory_pair_vector[0].first, 2);
    ASSERT_EQ(accessory_pair_vector[0].second, 2);

    ASSERT_EQ(accessory_pair_vector[1].first, 3);
    ASSERT_EQ(accessory_pair_vector[1].second, 10);

    ASSERT_EQ(core_pair_vector[0].first, 1);
    ASSERT_EQ(core_pair_vector[0].second, 14);
}

struct FrequencyExpectation {
    KmerType expected_kmer;
    std::vector<int> genome_idxs;
    std::vector<int> kmer_frequencies;
};

class KmerFrequencyClassifyTest : public testing::TestWithParam<FrequencyExpectation> {};

TEST_P(KmerFrequencyClassifyTest, testClassifyKmersForLoops) {
    const FrequencyExpectation& expectation = GetParam();
    auto kmer_classified = KmerUtils::classify_selected_kmers(expectation.kmer_frequencies, expectation.genome_idxs);
    ASSERT_EQ(expectation.expected_kmer, std::get<0>(kmer_classified));
}

INSTANTIATE_TEST_SUITE_P(testClassifyKmersForLoops, KmerFrequencyClassifyTest, testing::Values(
    FrequencyExpectation{
        KmerType::ACCESSORY,
        {1, 2, 3},
        {0, 14, 0, 11},
    },
    FrequencyExpectation{
        KmerType::CORE,
        {1, 2, 3},
        {0, 14, 2, 11},
    },
    FrequencyExpectation{
        KmerType::UNIQUE,
        {1, 2, 3},
        {1, 0, 0, 1},
    },
    FrequencyExpectation{
        KmerType::NOT_PRESENT,
        {1, 2, 3},
        {1, 0, 0, 0},
    }
));
TEST(testPangenomeStructure, testIterativelyCreatesLargerGenomeCombinations) {

}

TEST(testPangenomeStructure, testClassifyKmersForLoops) {

    int loops{5};
    int number_genomes{4};
    AnnotatedDBG* graph = nullptr;
    MockGraphWrapper wrapper = MockGraphWrapper(*graph);

    size_t frequencies{10000000};
    std::vector<std::vector<int>> freqs = create_test_frequencies_five_genomes(frequencies);

    boost::dynamic_bitset<> freqs_bitset = boost::dynamic_bitset<>(frequencies*5);
    boost::dynamic_bitset<>& bitset_ref{freqs_bitset};
    bitset_ref = fill_bitset_w_freqs(freqs, bitset_ref);

    vector<vector<int>> all_combinations = create_five_genome_test_combinations();

    std::vector<std::vector<int>> core_collection(5, std::vector<int>(3, 0));
    std::vector<std::vector<int>> accessory_collection(5, std::vector<int>(3, 0));
    std::vector<std::vector<int>> unique_collection(5, std::vector<int>(3, 0));
    PangenomeStructureAdapter::PangenomeStructureCollector collector{core_collection, accessory_collection, unique_collection, std::move(bitset_ref)};
    PangenomeStructureClassifier pangenome_classifier{wrapper, number_genomes, loops};

    collector = pangenome_classifier.loop_for_increasing_pangenome_size(collector, all_combinations);

    for (int loop = 0; loop < loops; loop++) {
        ASSERT_GT(std::reduce(collector.accessory_loops[loop].begin(), collector.accessory_loops[loop].end()) || std::reduce(collector.core_loops[loop].begin(), collector.core_loops[loop].end()) ||std::reduce(collector.unique_loops[loop].begin(), collector.unique_loops[loop].end()) , 0);
    }
}

TEST(testPangenomeStructure, testExportStructureMatrix) {
    size_t frequencies{100000};
    std::vector<std::vector<int>> freqs = create_test_frequencies_five_genomes(frequencies);
    vector<vector<int>> all_combinations = create_five_genome_test_combinations();
    boost::dynamic_bitset<> freqs_bitset = boost::dynamic_bitset<>(frequencies*5);
    boost::dynamic_bitset<>& bitset_ref{freqs_bitset};
    bitset_ref = fill_bitset_w_freqs(freqs, bitset_ref);

    int num_genomes{5};
    int num_loops{5};
    int combinations{num_genomes - 1};

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> distribution(0, 1000);

    std::vector<std::vector<int>> core_collection(num_loops, std::vector<int>(num_genomes, distribution(gen)));
    std::vector<std::vector<int>> accessory_collection(num_loops, std::vector<int>(num_genomes, distribution(gen)));
    std::vector<std::vector<int>> unique_collection(num_loops, std::vector<int>(num_genomes, distribution(gen)));

    PangenomeStructureAdapter::PangenomeStructureCollector collector{core_collection, accessory_collection, unique_collection, bitset_ref};
    auto kmer_matrix{*PangenomeStructureAdapter::export_matrix(collector)};
    for (int loop = 0; loop < num_loops; loop++) {
        for (int com_idx = 0; com_idx < combinations; com_idx++) {
            ASSERT_NE(kmer_matrix.core, nullptr);
            ASSERT_NE(kmer_matrix.accessory, nullptr);
            ASSERT_NE(kmer_matrix.unique, nullptr);
            ASSERT_EQ(kmer_matrix.core[combinations*loop + com_idx], core_collection[loop][com_idx]);
            ASSERT_EQ(kmer_matrix.accessory[combinations*loop + com_idx], accessory_collection[loop][com_idx]);
            ASSERT_EQ(kmer_matrix.unique[combinations*loop + com_idx], unique_collection[loop][com_idx]);
        }
    }
}

TEST(testPangenomeStructure, testDynamicBitSetforFreqs) {
    auto bs{boost::dynamic_bitset<>(10000000000)};
    bs.flip(0);
    ASSERT_EQ(bs[0], 1);
    std::cout << bs << std::endl;
    bs.flip(3);
    ASSERT_EQ(bs[3], 1);
    std::cout << bs << std::endl;
}