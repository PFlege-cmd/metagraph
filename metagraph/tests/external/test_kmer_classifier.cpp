//
// Created by Patrick Flege on 16/03/2026.
//

#include "../../../../../../../opt/homebrew/Cellar/boost/1.90.0_1/include/boost/iostreams/filter/zstd.hpp"
#include "cli/GraphWrapper.hpp"
#include "external/KmerClassifier.hpp"
#include "kmer/kmer.hpp"
#include "tests/cli/mock_graph_wrapper.h"
#include "tests/external/kmer_classifier_auxiliary.h"

#include "gmock/gmock-function-mocker.h"
#include "gtest/gtest.h"

using node_index_kmer = GraphWrapper::node_index_kmer;
using Label =  mtg::graph::AnnotatedSequenceGraph::Label;
using kmer_frequencies = GraphWrapper::kmer_frequencies;

TEST(testKmerClassifier, testGetGraph) {
    mtg::graph::AnnotatedDBG*  graph = nullptr;
    MockGraphWrapper wrapper = MockGraphWrapper(*graph);
    KmerClassifier kmer_classifier(wrapper);

    EXPECT_CALL(wrapper, get_number_nodes)
        .WillOnce(testing::Return(2));
    auto node_number = kmer_classifier.get_node_number();
    ASSERT_GT(node_number, 0);
}

TEST(testKmerClassifier, testGetCoordinates) {
    mtg::graph::AnnotatedDBG*  graph = nullptr;
    MockGraphWrapper wrapper = MockGraphWrapper(*graph);
    KmerClassifier kmer_classifier(wrapper);

    kmer_frequencies frequencies =  std::vector<std::tuple<Label, size_t, std::vector<size_t>>>();
    auto node_i = std::vector<node_index_kmer>();
    EXPECT_CALL(wrapper, get_kmer_frequencies(node_i)).WillOnce(testing::Return(frequencies));
    kmer_classifier.get_kmer_frequencies(node_i);
}

TEST(testKmerClassifier, testGetNumberOfGenomes) {
    mtg::graph::AnnotatedDBG*  graph = nullptr;
    MockGraphWrapper wrapper = MockGraphWrapper(*graph);

    int num_genomes_returned = 100;
    KmerClassifier kmer_classifier(wrapper);
    EXPECT_CALL(wrapper, get_num_genomes).WillOnce(testing::Return(num_genomes_returned));

    int num_genomes = kmer_classifier.get_num_genomes();
    ASSERT_GT(num_genomes, 0);
    ASSERT_EQ(num_genomes_returned, num_genomes);
}

TEST(testKmerClassifier, testGetGenomeNumberFromString) {
    mtg::graph::AnnotatedDBG*  graph = nullptr;
    MockGraphWrapper wrapper = MockGraphWrapper(*graph);

    KmerClassifier kmer_classifier(wrapper);
    std::string genome_contig_path("/test/my/contig_path/2_contig.fasta.gz");

    int genome_number = kmer_classifier.extractGenomeNumber(genome_contig_path);
    ASSERT_EQ(genome_number, 2);
}

TEST(testKmerClassifier, testGetGenomeNumberFromStringWithTripleDigit) {
    mtg::graph::AnnotatedDBG*  graph = nullptr;
    MockGraphWrapper wrapper = MockGraphWrapper(*graph);

    KmerClassifier kmer_classifier(wrapper);
    std::string genome_contig_path("/test/my/contig_path/202_contig.fasta.gz");

    int genome_number = kmer_classifier.extractGenomeNumber(genome_contig_path);
    ASSERT_EQ(genome_number, 202);
}

TEST(testKmerClassifier, testCalculateTotalCountMatrixOnce) {
    typedef std::tuple<Label, size_t, std::vector<size_t>> kmer_tuple;
    mtg::graph::AnnotatedDBG*  graph = nullptr;
    MockGraphWrapper wrapper = MockGraphWrapper(*graph);
    KmerClassifier kmer_classifier(wrapper, 3, 1, 4);

    Label genome_one_name = std::string("/0_contig/t");
    Label genome_two_name = std::string("/1_contig/t");
    Label genome_three_name = std::string("/2_contig");

    size_t genome_one_kmer_count = 3;
    size_t genome_two_kmer_count = 4;
    size_t genome_three_kmer_count = 5;

    auto first_genome_counts = std::make_tuple(genome_one_name, 1, std::vector<size_t>({genome_one_kmer_count}));
    auto second_genome_counts = std::make_tuple(genome_two_name, 1, std::vector<size_t>({genome_two_kmer_count}));
    auto third_genome_counts = std::make_tuple(genome_three_name, 1, std::vector<size_t>({genome_three_kmer_count}));

    kmer_frequencies frequencies_1 =  std::vector<kmer_tuple>({first_genome_counts, second_genome_counts, third_genome_counts});
    auto idx = std::vector<node_index_kmer>({1});

    EXPECT_CALL(wrapper, get_kmer_frequencies(idx)).WillOnce(testing::Return(frequencies_1));

    auto kmer_counts = kmer_classifier.count_kmer_per_genome(idx[0]);

    ASSERT_EQ(kmer_counts[0], 3);
    ASSERT_EQ(kmer_counts[1], 4);
    ASSERT_EQ(kmer_counts[2], 5);
}

TEST(testKmerClassifier, testThrowsExceptionAtInvalidGenomeNumber) {
    typedef std::tuple<Label, size_t, std::vector<size_t>> kmer_tuple;
    mtg::graph::AnnotatedDBG*  graph = nullptr;
    MockGraphWrapper wrapper = MockGraphWrapper(*graph);
    KmerClassifier kmer_classifier(wrapper);

    Label genome_one_name = std::string("/0_xxx/t");
    Label genome_two_name = std::string("/1_xxx/t");
    Label genome_three_name = std::string("/2_xxg");

    size_t genome_one_kmer_count = 3;
    size_t genome_two_kmer_count = 4;
    size_t genome_three_kmer_count = 5;

    auto first_genome_counts = std::make_tuple(genome_one_name, 1, std::vector<size_t>({genome_one_kmer_count}));
    auto second_genome_counts = std::make_tuple(genome_two_name, 1, std::vector<size_t>({genome_two_kmer_count}));
    auto third_genome_counts = std::make_tuple(genome_three_name, 1, std::vector<size_t>({genome_three_kmer_count}));

    kmer_frequencies frequencies_1 =  std::vector<kmer_tuple>({first_genome_counts, second_genome_counts, third_genome_counts});
    auto idx = std::vector<node_index_kmer>({1});

    EXPECT_CALL(wrapper, get_kmer_frequencies(idx)).WillOnce(testing::Return(frequencies_1));

    constexpr int size = 3;

    ASSERT_THROW(kmer_classifier.count_kmer_per_genome<size>(idx[0]), std::invalid_argument);
}

class TestKmerClassificationFixture : public ::testing::TestWithParam<std::tuple<KmerType, KmerClassifier::genomes_and_frequencies, std::array<int, 4>>> {};

INSTANTIATE_TEST_SUITE_P(
    IndividualKmerClassificationTests,
    TestKmerClassificationFixture,
    testing::Values(
        std::make_tuple<KmerType, KmerClassifier::genomes_and_frequencies, std::array<int, 4>>(
            KmerType::UNIQUE,
            std::vector<std::pair<int, int>>({
                std::make_pair(2, 5)}
            ),
            std::array<int, 4>({0, 0, 5, 0})
            ),
        std::make_tuple(
            KmerType::CORE,
            std::vector<std::pair<int, int>>(
            {
                std::make_pair(0, 1),
                std::make_pair(1, 2),
                std::make_pair(2, 1),
                std::make_pair(3, 2)
            }),
            std::array<int, 4>({1, 2, 1, 2})
        ),
        std::make_tuple(
            KmerType::ACCESSORY,
            std::vector<std::pair<int, int>>(
                {std::make_pair(0, 1),
                    std::make_pair(2, 1),
                    std::make_pair(3, 2)}
                ),
                std::array<int, 4>({1, 0, 1, 2})
                )
                ));


TEST_P(TestKmerClassificationFixture, testClassifiersKmersCorrectly) {
    auto values = GetParam();
    auto kmer_type = std::get<0>(values);
    auto occurrences= std::get<1>(values);
    auto kmercounts = std::get<2>(values);

    mtg::graph::AnnotatedDBG*  graph = nullptr;
    MockGraphWrapper wrapper = MockGraphWrapper(*graph);
    KmerClassifier kmer_classifier(wrapper, 4, 1);

    auto classification = kmer_classifier.classify_kmers(kmercounts);
    KmerType kmer_type_result = std::get<0>(classification);
    KmerClassifier::genomes_and_frequencies occurrences_result = std::get<1>(classification);

    ASSERT_EQ(kmer_type_result, kmer_type);
    ASSERT_EQ(occurrences_result, occurrences);

}

TEST(TestKmerClassification, testFillCountMatrix) {

    typedef std::tuple<Label, size_t, std::vector<size_t>> kmer_tuple;

    Label genome_one_name = std::string("/0_contig/t");
    Label genome_two_name = std::string("/1_contig/t");
    Label genome_three_name = std::string("/3_contig");

    size_t genome_one_kmer_count = 3;
    size_t genome_two_kmer_count = 4;
    size_t genome_three_kmer_count = 5;

    auto first_genome_counts = std::make_tuple(genome_one_name, 1, std::vector<size_t>({genome_one_kmer_count}));
    auto second_genome_counts = std::make_tuple(genome_two_name, 1, std::vector<size_t>({genome_two_kmer_count}));
    auto third_genome_counts = std::make_tuple(genome_three_name, 1, std::vector<size_t>({genome_three_kmer_count}));

    kmer_frequencies frequencies_1 =  std::vector<kmer_tuple>({first_genome_counts, second_genome_counts, third_genome_counts});

    Label genome_one_name_freq_2 = std::string("/0_contig/t");
    Label genome_two_name_freq_2 = std::string("/1_contig/t");
    Label genome_three_name_freq_2 = std::string("/2_contig");
    Label genome_four_name_freq_2 = std::string("/3_contig");

    size_t genome_one_kmer_count_freq_2 = 3;
    size_t genome_two_kmer_count_freq_2 = 4;
    size_t genome_three_kmer_count_freq_2 = 5;
    size_t genome_four_kmer_count_freq_2 = 5;

    auto first_genome_counts_freq_2 = std::make_tuple(genome_one_name_freq_2, 1, std::vector<size_t>({genome_one_kmer_count_freq_2}));
    auto second_genome_counts_freq_2 = std::make_tuple(genome_two_name_freq_2, 1, std::vector<size_t>({genome_two_kmer_count_freq_2}));
    auto third_genome_counts_freq_2 = std::make_tuple(genome_three_name_freq_2, 1, std::vector<size_t>({genome_three_kmer_count_freq_2}));
    auto fourth_genome_counts_freq_2 = std::make_tuple(genome_four_name_freq_2, 1, std::vector<size_t>({genome_four_kmer_count_freq_2}));

    kmer_frequencies frequencies_2 =  std::vector<kmer_tuple>({first_genome_counts_freq_2,
        second_genome_counts_freq_2,
        third_genome_counts_freq_2,
        fourth_genome_counts_freq_2});


    Label genome_two_name_freq_3 = std::string("/1_contig/t");

    size_t genome_two_kmer_count_freq_3 = 111;

    auto second_genome_counts_freq_3 = std::make_tuple(genome_two_name_freq_3, 1, std::vector<size_t>({genome_two_kmer_count_freq_3}));

    kmer_frequencies frequencies_3 =  std::vector<kmer_tuple>({
        second_genome_counts_freq_3});

    auto idx = std::vector<node_index_kmer>({1});
    auto idx2 = std::vector<node_index_kmer>({2});
    auto idx3 = std::vector<node_index_kmer>({3});
    auto idx4 = std::vector<node_index_kmer>({4});

    mtg::graph::AnnotatedDBG*  graph = nullptr;
    MockGraphWrapper wrapper = MockGraphWrapper(*graph);

    EXPECT_CALL(wrapper, get_number_nodes()).WillOnce(testing::Return(
        5));

    EXPECT_CALL(wrapper, get_kmer_frequencies(idx))
    .WillOnce(testing::Return(frequencies_1));
    EXPECT_CALL(wrapper, get_kmer_frequencies(idx2))
    .WillOnce(testing::Return(frequencies_2));
    EXPECT_CALL(wrapper, get_kmer_frequencies(idx3))
    .WillOnce(testing::Return(frequencies_1));
    EXPECT_CALL(wrapper, get_kmer_frequencies(idx4))
    .WillOnce(testing::Return(frequencies_3));

    KmerClassifier kmer_classifier(wrapper, 4, 1, 4);
    auto matrix = kmer_classifier.create_kmer_classification_matrix();

    ASSERT_EQ(matrix.size(), 3);
    ASSERT_EQ(matrix[1][3], 10);
    ASSERT_EQ(matrix[2][1], 111);
}

