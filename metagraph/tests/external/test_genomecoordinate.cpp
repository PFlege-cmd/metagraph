#include "external/GenomeCoordinate.h"
#include "tests/cli/mock_graph_wrapper.h"


#include "gtest/gtest.h"
//
// Created by Patrick Flege on 14/07/2026.
//
TEST(testGenomeCoordinateInitiation, testGenomeCoordinate) {
    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    std::string kmer = "ATGC";
    std::string genome = "genome1";
    auto sequence_lengths = std::vector<int>{10, 20, 30};
    GenomeCoordinate genome_coordinate(genome, kmer, sequence_lengths, wrapper);
    ASSERT_STREQ(genome_coordinate.get_reverse_kmer().c_str(), "GCAT");
}


TEST(testGenomeCoordinateLocate, testGenomeCoordinate) {
    SmallVector<uint64_t> occs = SmallVector<uint64_t>();
    occs.emplace_back(22);
    occs.emplace_back(33);
    auto coord_vector = std::vector<SmallVector<uint64_t>>({occs});
    GenomeCoordinate::genomeCoordinateTriple testCoordinateTriple = std::make_tuple<mtg::graph::AnnotatedDBG::Label, size_t, GenomeCoordinate::genomeCoordinateVector>("genome1", 10, std::move(coord_vector));

    std::vector<GenomeCoordinate::genomeCoordinateTriple> genomeCoordinateTriples = std::vector<GenomeCoordinate::genomeCoordinateTriple>({testCoordinateTriple});
    int first_coordinate_sequence_id = 2, second_coordinate_sequence_id = 3;
    int first_coordinate_position = 15, second_coordinate_position = 9;

    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    std::string kmer = "ATGC";
    std::string genome = "genome1";
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