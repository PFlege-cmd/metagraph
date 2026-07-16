#include <memory>

#include "mock_matrix_entry.h"
#include "asio/error.hpp"
#include "cli/DeBruijnGraphWrapper.h"
#include "cli/graph_glue.hpp"
#include "external/CoordinateRetriever.h"
#include "external/GenomeCoordinate.h"


#include "gtest/gtest.h"
#include "tests/cli/mock_graph_wrapper.h"

//
// Created by Patrick Flege on 14/07/2026.
//

using namespace std::string_view_literals;

struct HitsPerSequence {
    int* hitsPerSequence;
    int* sequenceOffsets;
    HitsPerSequence() = default;
    explicit HitsPerSequence(int* arr, int* offsets ) { hitsPerSequence = arr; sequenceOffsets = offsets; };
};
TEST(testRetrieveKmers, TestRegionOfInterest) {
    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    int max_frequency = 10;
    std::string genome = "genome1";
    int* seqs = new int[3] {10, 20, 50 };
    std::vector<std::unique_ptr<MatrixEntry>> genomeCoordinates;
    CoordinateRetriever* retriever
            = new CoordinateRetriever(3, max_frequency, genome, seqs,genomeCoordinates, wrapper);
    ASSERT_EQ(3, retriever->get_number_of_sequences());
}

TEST(testGetSequenceLength, TestRegionOfInterest) {
    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    int max_frequency = 10;
    std::string genome = "genome1";
    int* seqs = new int[3] {10, 20, 50 };
    std::vector<std::unique_ptr<MatrixEntry>> genomeCoordinates;

    CoordinateRetriever* retriever
            = new CoordinateRetriever(3, max_frequency, genome, seqs, genomeCoordinates, wrapper);

    ASSERT_EQ(10, retriever->get_sequence_length(0));
    ASSERT_EQ(20, retriever->get_sequence_length(1));
    ASSERT_EQ(50, retriever->get_sequence_length(2));
    ASSERT_THROW(retriever->get_sequence_length(-1), std::invalid_argument);
    ASSERT_THROW(retriever->get_sequence_length(3), std::invalid_argument);
}

TEST(testGetLocations, TestRegionOfInterest) {
    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    int max_frequency = 10;
    std::string genome = "genome1";
    int* seqs = new int[3] { 10, 20, 50 };
    std::string_view genome_view = std::string_view(genome);
    std::string_view kmer = std::string_view("ATGC");

    std::vector<std::array<int, 2>> coordinates_1 = std::vector({std::array<int, 2>{2, 11}, std::array<int, 2>{1, 11}});
    std::vector<std::array<int, 2>> coordinates_2 = std::vector({std::array<int, 2>{2, 11}, std::array<int, 2>{1, 11}});
    std::vector<std::array<int, 2>> coordinates_3 = std::vector({std::array<int, 2>{2, 11}, std::array<int, 2>{1, 11}});

    std::vector<int> lengts = std::vector<int>();

    MockMatrixEntry coord1 = MockMatrixEntry(genome_view, kmer, lengts, wrapper);
    MockMatrixEntry coord2 = MockMatrixEntry(genome_view, kmer, lengts, wrapper);
    MockMatrixEntry coord3 = MockMatrixEntry(genome_view, kmer, lengts, wrapper);

    std::unique_ptr<MockMatrixEntry> coord1_ptr  = std::make_unique<MockMatrixEntry>(coord1);
    std::unique_ptr<MockMatrixEntry> coord2_ptr  = std::make_unique<MockMatrixEntry>(coord2);
    std::unique_ptr<MockMatrixEntry> coord3_ptr  = std::make_unique<MockMatrixEntry>(coord3);

    auto c1 = coord1_ptr.get();
    auto c2 = coord2_ptr.get();
    auto c3 = coord3_ptr.get();

    std::vector<std::unique_ptr<MatrixEntry>> coordinates = std::vector<std::unique_ptr<MatrixEntry>>();
    coordinates.push_back(std::move(coord1_ptr));
    coordinates.push_back(std::move(coord2_ptr));
    coordinates.push_back(std::move(coord3_ptr));


    //EXPECT_CALL(wrapper, get_K()).WillRepeatedly(testing::Return(kmer_size));
    EXPECT_CALL(*c1, locate_both_kmers).WillOnce(testing::Return(coordinates_1));
    EXPECT_CALL(*c2, locate_both_kmers).WillOnce(testing::Return(coordinates_2));
    EXPECT_CALL(*c3, locate_both_kmers).WillOnce(testing::Return(coordinates_3));
    CoordinateRetriever* retriever
            = new CoordinateRetriever(3, max_frequency, genome, seqs, coordinates, wrapper);

    //retriever->setCoordinates(coordinates);
    std::vector<std::array<int, 2>> kmer_coords = retriever->get_kmer_positions();
}

TEST(testIntegration, TestRegionOfInterest) {
    std::string data_path = std::string("../../../patrick-pan-tools/chloroplast_DB");

    graph_glue glue = graph_glue(0, NULL);

    //static std::shared_ptr<AnnotatedDBG> graph_coord = glue.load_coord_dbg(data_path);


    std::string_view kmer1("ATCAACCGAGGTA");
    std::string_view kmer2("AAGCAAGGAAAAA");
    std::string_view kmer3("AATCAACCCGAAG");
    std::string_view kmer4("GCAGGGTCTTCCA");
    std::string_view kmer5("AGCATTCAATTAG");
    std::string_view kmer6("CTTGGGAAGATCA");
    std::string_view kmer8("AATCCTAAAAAAA");
    std::string_view kmer9("CTAATTCGGCAAA");
    std::string_view kmer10("AATTGTCAGCAAA");
    std::string_view kmer11("AATACAACACCAA");
    std::string_view kmer12("ATAAAACGATAAA");
    std::string_view kmer13("CTAAATTTATTGA");
    std::string_view kmer14("ACCCAATTTGGGA");
    std::string_view kmer15("AAGATATTGGAAC");
    std::string_view kmer16("ACAAAAAGGGTTA");
    std::string_view kmer17("AATATTTTATTAA");
    std::string_view kmer18("TAAAATTTCCAGA");
    std::string_view kmer19("AAATGGTAGCAAC");
    std::string_view kmer20("AATGGATAGGACA");

    std::vector<std::string_view> kmers = {kmer1, kmer2, kmer3, kmer4, kmer5, kmer6, kmer8, kmer9, kmer10, kmer11, kmer12, kmer13, kmer14, kmer15, kmer16, kmer17, kmer18, kmer19, kmer20};
    static std::shared_ptr<AnnotatedDBG> graph = glue.load_coord_dbg(data_path);

    auto wrapper = DeBruijnGraphWrapper(*graph);
    //std::string_view kmer = std::string_view("ACCCAATTTGGGA");
    std::string_view genome_view = std::string_view("chloroplasts/O_sativa.fasta");
    std::vector<int> lengths = std::vector<int>({ 134496 });

    std::vector<unique_ptr<MatrixEntry>> coordinates
            = std::vector<std::unique_ptr<MatrixEntry>>();
    for_each(kmers.begin(), kmers.end(), [&coordinates, &genome_view, &wrapper, &lengths](auto& kmer) {
        coordinates.push_back(std::make_unique<GenomeCoordinate>(genome_view, kmer, lengths, wrapper));
    });

    //int genome_number = 3;
    int* seqs = new int[3] { 10, 20, 50 };
    std::string genome = std::string("genome");
    //auto genomeCoordinate = new GenomeCoordinate(genome_view, kmer, lengths, wrapper);
    //std::vector<std::array<int, 2>> results = genomeCoordinate->locate_both_kmers();
    CoordinateRetriever* retriever = new CoordinateRetriever(1, 5, genome, seqs, coordinates, wrapper);
    std::vector<std::array<int, 2>> results = retriever->get_kmer_positions();

    ASSERT_EQ(3, results.size());

}



