//
// Created by Patrick Flege on 30/06/2026.
//

#include "cli/DeBruijnGraphWrapper.h"
#include "cli/graph_glue.hpp"
#include "external/ReadMapper.h"
#include "tests/cli/mock_graph_wrapper.h"

#include "gtest/gtest.h"
struct HitsPerSequence {
    int* hitsPerSequence;
    int* sequenceOffsets;
    HitsPerSequence() = default;
    explicit HitsPerSequence(int* arr, int* offsets ) { hitsPerSequence = arr; sequenceOffsets = offsets; };
};
TEST(testReadmapping, testInvalidRead) {
    std::string db_string = std::string("/Users/patrick_flege/git/patrick-pan-tools/chloroplast_DB/");
    std::string gd_string =  std::string("chloroplasts/O_sativa.fasta");
    char * database  = db_string.data();
    char * genome_name = gd_string.data();

    std::array<int, 2> sq_lengths_array = std::array<int, 2>{ 0, 155293 };
    int * sequence_lengths = sq_lengths_array.data();
    int no_of_sequences = 2;

    std::string readchar = std::string("AAGCCTATGGGGTCGTTTCTGTACCTGGATAACCAGCACTGAGAACCGTCTTTACATTGGATGGTTTGGT");
    std::string readchar_fourth_line_c_sativus{"AAAGGAGCAATACCAATCCTCTTGAGAAAACAAGAAATTGGTTATTGCTCCTTTACTTTCAAGAACTCGT"}; // As name says, fourth line in C sativus fasta
    char * read = readchar.data();
    graph_glue glue = graph_glue(0, NULL);
    static std::shared_ptr<AnnotatedDBG> graph = glue.load_coord_dbg(database);
    HitsPerSequence* results = new HitsPerSequence();
    glue.do_pantools_work(genome_name, sequence_lengths, no_of_sequences, results, read, graph);
    //return results;
}


TEST(testReadmapping, testInvalidReadClassBased) {

    int* seqsPerGenome = new int[4]{3, 1, 1, 1};
    int* seqLengths = new int[6]{ 155293,155293, 155293, 134496, 155461, 155296};
    std::string db_string = std::string("/Users/patrick_flege/git/patrick-pan-tools/chloroplast_DB_c_sativus_tripled/");
    std::string gd_string =  std::string("chloroplasts/O_sativa.fasta");
    char * database  = db_string.data();
    std::string read = std::string("AAGCCTATGGGGTCGTTTCTGTACCTGGATAACCAGCACTGAGAACCGTCTTTACATTGGATGGTTTGGT");
    std::string readchar_fourth_line_c_sativus{"AAAGGAGCAATACCAATCCTCTTGAGAAAACAAGAAATTGGTTATTGCTCCTTTACTTTCAAGAACTCGT"}; // As name says, fourth line in C sativus fasta

    //using sequence_and_position = ReadMapper::sequence_and_position;

        graph_glue glue = graph_glue(0, NULL);
        static std::shared_ptr<AnnotatedDBG> graph = glue.load_coord_dbg(database);
        auto wrapper = DeBruijnGraphWrapper(*graph);


        ReadMapper read_mapper(wrapper);

        read_mapper.initialize_all_fields(seqLengths, seqsPerGenome);
        read_mapper.set_read(std::string_view{readchar_fourth_line_c_sativus});

        read_mapper.read_mapping_forward();
        read_mapper.read_mapping_backward();

        NativeStructures::HitsPerGenomes* hitsPerGenome {read_mapper.transform_to_native_output()};
        ASSERT_NE(nullptr, hitsPerGenome);
        ASSERT_GT(hitsPerGenome->hitFreqPerGenome[0], 10);
        ASSERT_EQ(hitsPerGenome->hitFreqPerGenome[1], 0);
        ASSERT_GT(hitsPerGenome->hitsPerSequence[0].sequenceOffsets[0], 0);
        ASSERT_GT(hitsPerGenome->hitsPerSequence[0].sequenceOffsets[1], 0);
        ASSERT_GT(hitsPerGenome->hitsPerSequence[0].sequenceOffsets[2], 0);
}

TEST(testReadmapping, testConvertToProperOutput) {
    std::vector<std::vector<unsigned long long>> sequenceOffsets{};
    std::vector<std::vector<unsigned long long>> sequenceLengths{{100, 102}, {120, 2000}, {1000, 2000, 3000, 50000, 50000}};
    std::map<std::string_view, std::vector<std::array<int, 2>>> genome_map{};
    std::map<std::string_view, size_t> genome_name_to_index_map{};
    mtg::graph::AnnotatedDBG* graph{nullptr};

    MockGraphWrapper wrapper{*graph};

    std::string read{"AGGGGGGGAGGGGGGA"};

    genome_map.emplace("Agull", std::vector<std::array<int, 2>>{{1, 1}, {1, 20}, {2, 200}});
    genome_map.emplace("Smeagol", std::vector<std::array<int, 2>>{{1, 2}, {1, 30}, {2, 1000}});
    genome_map.emplace("Kazadhum", std::vector<std::array<int, 2>>{{1, 3}, {4, 42}, {5, 6666}});

    genome_name_to_index_map.emplace("Agull", 0);
    genome_name_to_index_map.emplace("Smeagol", 1);
    genome_name_to_index_map.emplace("Kazadhum", 2);

    ReadMapper readMapper(
            sequenceOffsets,
            sequenceLengths,
            genome_map,
            genome_name_to_index_map,
            wrapper,
            std::string_view{read}
            );

    NativeStructures::HitsPerGenomes* hitsPerGenome {readMapper.transform_to_native_output()};

    NativeStructures::HitsPerSequence h1 = hitsPerGenome->hitsPerSequence[0];
    NativeStructures::HitsPerSequence* h2 = &hitsPerGenome->hitsPerSequence[1];
    NativeStructures::HitsPerSequence* h3 = &hitsPerGenome->hitsPerSequence[2];

    // use genome_to_index map
    ASSERT_EQ(h1.sequenceOffsets[0], 2);
    ASSERT_EQ(h1.sequenceOffsets[1], 1);

    ASSERT_EQ(h2->sequenceOffsets[0], 2);
    ASSERT_EQ(h2->sequenceOffsets[1], 1);

    ASSERT_EQ(h2->hitsPerSequence[0], 2);
    ASSERT_EQ(h2->hitsPerSequence[1], 30);
    ASSERT_EQ(h2->hitsPerSequence[2], 1000);

    ASSERT_EQ(h3->sequenceOffsets[0], 1);
    ASSERT_EQ(h3->sequenceOffsets[1], 0);
    ASSERT_EQ(h3->sequenceOffsets[2], 0);
    ASSERT_EQ(h3->sequenceOffsets[3], 1);
    ASSERT_EQ(h3->sequenceOffsets[4], 1);

    ASSERT_EQ(h3->hitsPerSequence[0], 3);
    ASSERT_EQ(h3->hitsPerSequence[1], 42);
    ASSERT_EQ(h3->hitsPerSequence[2], 6666);

    ASSERT_EQ(hitsPerGenome->hitFreqPerGenome[0], 3);
    delete hitsPerGenome;
}