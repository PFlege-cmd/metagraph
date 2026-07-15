//
// Created by Patrick Flege on 30/06/2026.
//

#include "cli/graph_glue.hpp"

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
    char * read = readchar.data();
    graph_glue glue = graph_glue(0, NULL);
    static std::shared_ptr<AnnotatedDBG> graph = glue.load_coord_dbg(database);
    HitsPerSequence* results = new HitsPerSequence();
    glue.do_pantools_work(genome_name, sequence_lengths, no_of_sequences, results, read, graph);
    //return results;
}