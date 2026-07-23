#include <memory>

#include "mock_matrix_entry.h"
#include "mock_matrix_entry_factory.h"
#include "asio/error.hpp"
#include "cli/DeBruijnGraphWrapper.h"
#include "cli/graph_glue.hpp"
#include "external/CoordinateRetriever.h"
#include "external/GenomeCoordinate.h"
#include "external/GenomeCoordinateFactory.h"
#include "external/KmerClassifier.hpp"


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
            = new CoordinateRetriever(3, max_frequency, genome, seqs, wrapper);
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
            = new CoordinateRetriever(3, max_frequency, genome, seqs, wrapper);

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
            = new CoordinateRetriever(3, max_frequency, genome, seqs, wrapper);

    //retriever->setCoordinates(coordinates);
    std::vector<std::array<int, 2>> kmer_coords = retriever->get_kmer_positions(coordinates);
}

TEST(testCreatesEntries, TestRegionOfInterest) {
    mtg::graph::AnnotatedDBG* graph = nullptr;
    auto wrapper = MockGraphWrapper(*graph);
    int max_frequency = 10;
    std::string genome = "genome1";
    int* seqs = new int[3] { 10, 20, 50 };
    int no_seqs = 3;
    std::string_view genome_view = std::string_view(genome);
    char const* one = "ATGC";
    char const* two = "ATGC";
    char const* three = "ATGC";
    char const* four = "ATGC";
    char const* five = "ATGC";

    char const* kmer_raw[5] = {one, two, three, four, five};
    std::string_view kmer = std::string_view("ATGC");
    std::vector<int> lengts = std::vector<int>();

    // std::string_view& genome,
    //                               std::string_view& kmer,
    //                               std::vector<int>& sequence_lengths,
    //                               GraphWrapper& graph
    //auto matrixEntryFactory = new MockMatrixEntryFactory();
    //MatrixEntryFactory::instance = matrixEntryFactory;

    auto ptr1 = std::make_unique<MockMatrixEntry>(genome_view, kmer,lengts, wrapper);
    auto ptr2 = std::make_unique<MockMatrixEntry>(genome_view, kmer,lengts, wrapper);
    auto ptr3 = std::make_unique<MockMatrixEntry>(genome_view, kmer,lengts, wrapper);
    std::vector<std::unique_ptr<MatrixEntry>> coordinates = std::vector<std::unique_ptr<MatrixEntry>>();
    coordinates.push_back(std::move(ptr1));
    coordinates.push_back(std::move(ptr2));
    coordinates.push_back(std::move(ptr3)); // If not moved, will not work, since this is a unique pointer.
    auto inst =
            MockMatrixEntryFactory::getInstance(genome, no_seqs, seqs, wrapper).get(); // Ohhh... WIll attempt copy-elision

    EXPECT_CALL(*inst, createKmerEntries).WillOnce(testing::Return(testing::ByMove(std::move(coordinates))));

    CoordinateRetriever* retriever = new CoordinateRetriever(3, max_frequency, genome, seqs, wrapper);
    std::vector<std::unique_ptr<MatrixEntry>> entries = retriever->createEntries(*inst, kmer_raw, 5);
    ASSERT_GE(3, entries.size());
    delete inst;
}

TEST(testFilterKmers, TestRegionOfInterest) {
    typedef std::tuple<std::string, size_t, std::vector<size_t>> kmer_tuple;
    typedef std::string Label;
    mtg::graph::AnnotatedDBG* graph = nullptr;
    int max_frequency = 5;
    size_t expected_kmers = 4;
    std::string genome = "genome1";

    std::string_view kmer1("ATCAACCGAGGTA");
    std::string_view kmer2("AAGCAAGGAAAAA");
    std::string_view kmer3("AATCAACCCGAAG");
    std::string_view kmer4("GCAGGGTCTTCCA");
    std::string_view kmer5("GCAGGGTATTCCA");

    std::vector<string_view> kmers = {kmer1, kmer2, kmer3, kmer4, kmer5};

    Label genome_one_name = std::string("/0_contig/t");
    Label genome_two_name = std::string("/1_contig/t");
    Label genome_three_name = std::string("/2_contig");

    size_t genome_one_kmer_count = 3;
    size_t genome_one_kmer_count_high = 3;
    size_t genome_two_kmer_count = 1;
    size_t genome_two_kmer_count_high = 4;

    size_t genome_three_kmer_count = 1;
    size_t genome_three_kmer_count_high = 5;
    size_t genome_three_kmer_count_too_much = 6;

    auto first_genome_counts = std::make_tuple(genome_one_name, 1, std::vector<size_t>({genome_one_kmer_count}));
    auto first_genome_counts_high = std::make_tuple(genome_one_name, 1, std::vector<size_t>({genome_one_kmer_count_high}));

    auto second_genome_counts = std::make_tuple(genome_two_name, 1, std::vector<size_t>({genome_two_kmer_count}));
    auto second_genome_counts_high = std::make_tuple(genome_two_name, 1, std::vector<size_t>({genome_two_kmer_count_high}));

    auto third_genome_counts = std::make_tuple(genome_three_name, 1, std::vector<size_t>({genome_three_kmer_count}));
    auto third_genome_counts_high = std::make_tuple(genome_three_name, 1, std::vector<size_t>({genome_three_kmer_count_high}));
    auto third_genome_counts_too_much = std::make_tuple(genome_three_name, 1, std::vector<size_t>({genome_three_kmer_count_too_much}));


    kmer_frequencies frequencies_1 =  std::vector<kmer_tuple>({first_genome_counts, second_genome_counts, third_genome_counts});
    kmer_frequencies frequencies_2 =  std::vector<kmer_tuple>({first_genome_counts, second_genome_counts, third_genome_counts_too_much});
    kmer_frequencies frequencies_3 =  std::vector<kmer_tuple>({first_genome_counts_high, second_genome_counts_high, third_genome_counts_high});


    int* seqs = new int[3] { 10, 20, 50 };
    auto wrapper = MockGraphWrapper(*graph);
    auto wrapper2 = MockGraphWrapper(*graph);
    EXPECT_CALL(wrapper, get_kmer_frequencies(kmer1)).WillOnce(testing::Return(frequencies_1));
    EXPECT_CALL(wrapper, get_kmer_frequencies(kmer2)).WillOnce(testing::Return(frequencies_1));
    EXPECT_CALL(wrapper, get_kmer_frequencies(kmer3)).WillOnce(testing::Return(frequencies_1));
    EXPECT_CALL(wrapper, get_kmer_frequencies(kmer4)).WillOnce(testing::Return(frequencies_2));
    EXPECT_CALL(wrapper, get_kmer_frequencies(kmer5)).WillOnce(testing::Return(frequencies_3));


    auto retriever = CoordinateRetriever(4, max_frequency, genome_three_name, seqs, wrapper);
    auto max_filtered = retriever.filter_max_frequency(kmers);
    ASSERT_EQ(max_filtered.size(), expected_kmers);
}

TEST(testIntegration, TestRegionOfInterest) {
    std::string data_path = std::string("../../../patrick-pan-tools/chloroplast_DB");

    graph_glue glue = graph_glue(0, NULL);
    int genome_nr = 0;

    //static std::shared_ptr<AnnotatedDBG> graph_coord = glue.load_coord_dbg(data_path);
    //int * seqs_ptr = new int[1]{134496};
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

    const char * raw_kmers[] = {"ATCAACCGAGGTA", "AAGCAAGGAAAAA", "AATCAACCCGAAG", "GCAGGGTCTTCCA"};
    int no_kmers = sizeof(raw_kmers) / sizeof(raw_kmers[0]);
    const char * testykmer = "act";
    std::cout << raw_kmers[2] << std::endl;
    std::cout << testykmer << std::endl;
    std::cout << no_kmers << std::endl;
            std::vector<std::string_view> kmers = {kmer1, kmer2, kmer3, kmer4, kmer5, kmer6, kmer8, kmer9, kmer10, kmer11, kmer12, kmer13, kmer14, kmer15, kmer16, kmer17, kmer18, kmer19, kmer20};
    static std::shared_ptr<AnnotatedDBG> graph = glue.load_coord_dbg(data_path);
    static std::shared_ptr<AnnotatedDBG> kmer_graph = glue.load_dbg(data_path);

    auto wrapper = DeBruijnGraphWrapper(*graph);
    auto kmer_graph_wrapper = DeBruijnGraphWrapper(*kmer_graph);
    //std::string_view kmer = std::string_view("ACCCAATTTGGGA");

    auto labels = kmer_graph_wrapper.get_graph()->get_annotator().get_label_encoder().get_labels();

    auto kmer_labels = kmer_graph_wrapper.get_graph()->get_annotator().get_label_encoder().get_labels();
    auto coord_labels = wrapper.get_graph()->get_annotator().get_label_encoder().get_labels();

    std::string genome_string = kmer_labels[genome_nr];
    std::string_view genome_view_string = std::string_view(coord_labels[genome_nr]);

    //std::string_view genome_view = std::string_view("chloroplasts/O_sativa.fasta");
    std::vector<int> lengths = std::vector<int>({ 134496 });

    //auto factory = GenomeCoordinateFactory::getInstance("chloroplasts/O_sativa.fasta", 2, seqs_ptr, wrapper).get();
    //factory->createKmerEntries(raw_kmers, no_kmers);
    int* seqs = new int[3] { 10, 20, 50 };
    //std::string genome = std::string("chloroplasts/O_sativa.fasta");
    CoordinateRetriever* retriever = new CoordinateRetriever(1, 5, genome_string, seqs, kmer_graph_wrapper);
    kmers = retriever->filter_max_frequency(kmers);

    retriever->set_graph(wrapper);

    std::vector<unique_ptr<MatrixEntry>> coordinates
            = std::vector<std::unique_ptr<MatrixEntry>>();
    for_each(kmers.begin(), kmers.end(), [&coordinates, &genome_view_string, &wrapper, &lengths](auto& kmer) {
        coordinates.push_back(std::make_unique<GenomeCoordinate>(genome_view_string, kmer, lengths, wrapper));
    });


    std::vector<std::array<int, 2>> results = retriever->get_kmer_positions(coordinates);

    ASSERT_EQ(3, results.size());
}



