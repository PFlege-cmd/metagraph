#include "external/GenomeCoordinate.h"
#include "tests/cli/mock_graph_wrapper.h"


#include "gtest/gtest.h"
using namespace std::string_view_literals;
//
// Created by Patrick Flege on 14/07/2026.
//
TEST(testGenomeCoordinateInitiation, testGenomeCoordinate) {
    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    std::string_view kmer = "ATGC";
    std::string_view genome = "genome1";
    auto sequence_lengths = std::vector<int>{10, 20, 30};
    GenomeCoordinate genome_coordinate(genome, kmer, sequence_lengths, wrapper);
    ASSERT_STREQ(genome_coordinate.get_reverse_kmer().data(), "GCAT");
}


TEST(testGenomeCoordinateLocate, testGenomeCoordinate) {
    SmallVector<uint64_t> occs = SmallVector<uint64_t>();
    occs.emplace_back(22);
    occs.emplace_back(33);
    auto coord_vector = std::vector<SmallVector<uint64_t>>({occs});
    GenomeCoordinate::genomeCoordinateTriple testCoordinateTriple = std::make_tuple<mtg::graph::AnnotatedDBG::Label, size_t, GenomeCoordinate::genomeCoordinateVector>("genome1", 10, std::move(coord_vector));

    std::vector<GenomeCoordinate::genomeCoordinateTriple> genomeCoordinateTriples = std::vector<GenomeCoordinate::genomeCoordinateTriple>({testCoordinateTriple});
    int first_coordinate_sequence_id = 2, second_coordinate_sequence_id = 3;
    int first_coordinate_position = 16, second_coordinate_position = 10;

    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    std::string_view kmer = "ATGC";
    std::string_view genome = "genome1";
    auto sequence_lengths = std::vector<int>{10, 20, 30};
    size_t kmer_size = 4;
    EXPECT_CALL(wrapper, get_kmer_coordinates(kmer))
        .WillOnce(testing::Return(genomeCoordinateTriples));
    EXPECT_CALL(wrapper, get_K()).WillRepeatedly(testing::Return(kmer_size));

    GenomeCoordinate genome_coordinate(genome, kmer, sequence_lengths, wrapper);
    std::vector<std::array<int, 2>> coordinate = genome_coordinate.locate();
    std::array<int, 2> coord_1 = coordinate[0], coord_2 = coordinate[1];

    ASSERT_EQ(coord_1[0], first_coordinate_sequence_id);
    ASSERT_EQ(coord_2[0], second_coordinate_sequence_id);
    ASSERT_EQ(coord_1[1], first_coordinate_position);
    ASSERT_EQ(coord_2[1], second_coordinate_position);
}

TEST(testGenomeCoordinateLocateReverse, testGenomeCoordinate) {
    SmallVector<uint64_t> occs = SmallVector<uint64_t>();
    occs.emplace_back(22);
    occs.emplace_back(33);
    auto coord_vector = std::vector<SmallVector<uint64_t>>({occs});
    GenomeCoordinate::genomeCoordinateTriple testCoordinateTriple = std::make_tuple<mtg::graph::AnnotatedDBG::Label, size_t, GenomeCoordinate::genomeCoordinateVector>("genome1", 10, std::move(coord_vector));

    std::vector<GenomeCoordinate::genomeCoordinateTriple> genomeCoordinateTriples = std::vector<GenomeCoordinate::genomeCoordinateTriple>({testCoordinateTriple});
    int first_coordinate_sequence_id = 2, second_coordinate_sequence_id = 3;
    int first_coordinate_position = 20, second_coordinate_position = 14;

    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    std::string_view kmer = "ATGC";
    std::string_view rev_kmer = "GCAT";
    std::string_view genome = "genome1";
    auto sequence_lengths = std::vector<int>{10, 20, 30};
    size_t kmer_size = 4;
    EXPECT_CALL(wrapper, get_kmer_coordinates(rev_kmer))
        .WillOnce(testing::Return(genomeCoordinateTriples));
    EXPECT_CALL(wrapper, get_K()).WillRepeatedly(testing::Return(kmer_size));

    GenomeCoordinate genome_coordinate(genome, kmer, sequence_lengths, wrapper);
    std::vector<std::array<int, 2>> coordinate = genome_coordinate.locate_reverse();
    std::array<int, 2> coord_1 = coordinate[0], coord_2 = coordinate[1];

    ASSERT_EQ(coord_1[0], first_coordinate_sequence_id);
    ASSERT_EQ(coord_2[0], second_coordinate_sequence_id);
    ASSERT_EQ(coord_1[1], first_coordinate_position);
    ASSERT_EQ(coord_2[1], second_coordinate_position);
}

TEST(testGenomeCoordinateEmpty, testGenomeCoordinate) {
    SmallVector<uint64_t> occs = SmallVector<uint64_t>();
    occs.emplace_back(22);
    occs.emplace_back(33);
    auto coord_vector = std::vector<SmallVector<uint64_t>>({occs});
    GenomeCoordinate::genomeCoordinateTriple testCoordinateTriple = std::make_tuple<mtg::graph::AnnotatedDBG::Label, size_t, GenomeCoordinate::genomeCoordinateVector>("genome2", 10, std::move(coord_vector));

    std::vector<GenomeCoordinate::genomeCoordinateTriple> genomeCoordinateTriples = std::vector<GenomeCoordinate::genomeCoordinateTriple>({testCoordinateTriple});

    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    std::string_view kmer = "ATGC";
    std::string_view rev_kmer = "GCAT";
    std::string_view genome = "genome1";
    auto sequence_lengths = std::vector<int>{10, 20, 30};
    size_t kmer_size = 4;
    EXPECT_CALL(wrapper, get_kmer_coordinates(rev_kmer))
        .WillOnce(testing::Return(genomeCoordinateTriples));
    EXPECT_CALL(wrapper, get_K()).WillRepeatedly(testing::Return(kmer_size));

    GenomeCoordinate genome_coordinate(genome, kmer, sequence_lengths, wrapper);
    std::vector<std::array<int, 2>> coordinate = genome_coordinate.locate_reverse();

    ASSERT_EQ(0, coordinate.size());
}

TEST(testFindsAllMatches, testGenomeCoordinate) {
    SmallVector<uint64_t> occs = SmallVector<uint64_t>();
    occs.emplace_back(22);
    occs.emplace_back(33);
    auto coord_vector = std::vector<SmallVector<uint64_t>>({occs});
    auto coord_vector_2 = std::vector<SmallVector<uint64_t>>({occs});
    GenomeCoordinate::genomeCoordinateTriple testCoordinateTriple = std::make_tuple<mtg::graph::AnnotatedDBG::Label, size_t, GenomeCoordinate::genomeCoordinateVector>("genome1", 10, std::move(coord_vector));
    GenomeCoordinate::genomeCoordinateTriple testCoordinateTripleReverse = std::make_tuple<mtg::graph::AnnotatedDBG::Label, size_t, GenomeCoordinate::genomeCoordinateVector>("genome1", 10, std::move(coord_vector_2));

    std::vector<GenomeCoordinate::genomeCoordinateTriple> genomeCoordinateTriples = std::vector<GenomeCoordinate::genomeCoordinateTriple>({testCoordinateTriple});
    std::vector<GenomeCoordinate::genomeCoordinateTriple> genomeCoordinateTriplesReverse = std::vector<GenomeCoordinate::genomeCoordinateTriple>({testCoordinateTripleReverse});

    int first_coordinate_sequence_id = 2, second_coordinate_sequence_id = 3;
    int first_coordinate_position = 16, second_coordinate_position = 10;

    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    auto kmer = "ATGC"sv;
    auto rev_kmer = "GCAT"sv;
    auto genome = "genome1"sv;
    auto sequence_lengths = std::vector<int>{10, 20, 30};
    size_t kmer_size = 4;
    EXPECT_CALL(wrapper, get_kmer_coordinates(kmer))
        .WillOnce(testing::Return(genomeCoordinateTriples));
    EXPECT_CALL(wrapper, get_kmer_coordinates(rev_kmer))
        .WillOnce(testing::Return(genomeCoordinateTriplesReverse));

    EXPECT_CALL(wrapper, get_K()).WillRepeatedly(testing::Return(kmer_size));

    GenomeCoordinate genome_coordinate(genome, kmer, sequence_lengths, wrapper);
    std::vector<std::array<int, 2>> coordinate = genome_coordinate.locate_both_kmers();
    std::array<int, 2> coord_1 = coordinate[0], coord_2 = coordinate[1], coord_3 = coordinate[2], coord_4 = coordinate[3];

    ASSERT_EQ(4, coordinate.size());
    ASSERT_EQ(first_coordinate_sequence_id, coord_1[0]);
    ASSERT_EQ(second_coordinate_sequence_id, coord_2[0]);
    ASSERT_EQ(first_coordinate_position, coord_1[1]);
    ASSERT_EQ(second_coordinate_position, coord_2[1]);

    ASSERT_EQ(first_coordinate_sequence_id, coord_3[0]);
    ASSERT_EQ(second_coordinate_sequence_id, coord_4[0]);
    ASSERT_EQ(first_coordinate_position + kmer_size, coord_3[1]) ;
    ASSERT_EQ(second_coordinate_position + kmer_size, coord_4[1]);
}