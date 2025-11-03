//
// Created by Patrick Flege on 23/07/2025.
//

#include "load/load_annotated_graph.hpp"
#include <filesystem>
#include <iostream>

#include "load/load_graph.hpp"
#include "graph_glue.hpp"

#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
using namespace mtg::graph;
namespace fs = std::filesystem;

extern "C"{
    __attribute__((visibility("default")))
    __attribute__((used))
    int get_very_special_secret_message() {
        load_dbg();
        return 42;
    }
}

extern "C"{
    __attribute__((visibility("default")))
    __attribute__((used))
    int receive_and_return_double_secret_message(const int x) {
        int y = 2*x;
        return y;
    }
}
struct HitsPerSequence {
    int* hitsPerSequence;
    int* sequenceOffsets;
    HitsPerSequence() = default;
    explicit HitsPerSequence(int* arr, int* offsets ) { hitsPerSequence = arr; sequenceOffsets = offsets; };
};



extern "C"{
    HitsPerSequence* retrieve_hits_for_genome(char * genome_name, int* sequence_lengths, int no_of_sequences, char * read) {
        for (int i = 0; i < no_of_sequences; i++) {
        std::cout << sequence_lengths[i] << std::endl;
    }

    HitsPerSequence* results = new HitsPerSequence();
    do_pantools_work(genome_name, sequence_lengths, no_of_sequences, results, read);

    return results;

    }

    void free_hits(HitsPerSequence* hits, int size_of_array) {
        for (int i = 0; i < size_of_array; i++) {
            delete hits[i].hitsPerSequence;
            delete hits[i].sequenceOffsets;
        }
        delete[] hits;
    }
}





extern "C"{
    __attribute__((visibility("default")))
    __attribute__((used))
    const char *  retrieve_sequence_for_coordinates(char * genome_name, int start, int end) {
        static std::shared_ptr<AnnotatedDBG> graph = load_dbg();
        std::string seq = graph->get_sequence_for_coords(genome_name, start, end);
        std::cout << "Checking genome: " << std::endl;
        std::cout << genome_name << std::endl;
        char * result = (char*)malloc(seq.length() + 1);
        std::strcpy(result, seq.c_str());
        return result;
    }
}

extern "C"{
    const char *  retrieve_kmer_for_coordinates(int genome, int start) {
        std::cout << "Checking genome: " << std::endl;
        std::cout << genome << std::endl;
        std::cout << start << std::endl;
        static std::shared_ptr<AnnotatedDBG> graph = load_dbg();
        std::string genome_name = graph->get_annotator().get_label_encoder().get_labels()[genome];
        std::string kmer = graph->get_kmer_for_coords(genome_name, start);
        char * result = (char*)malloc(kmer.length() + 1);
        std::strcpy(result, kmer.c_str());
        return result;
    }

    void free_memory(const char* str) {
        free((void*)str);
    }
}

void do_pantools_work(char* genome_name,
             int* sequence_lengths,
             int no_of_sequences,
             HitsPerSequence*& results, char* read_ptr) {
    static std::shared_ptr<AnnotatedDBG> graph = load_dbg();
    std::string_view read = read_ptr;
        //"AGTACCAGAGATTCCTAGAGGCATAC";
    // std::string_view genome = "../succinct_data/zapdos.fasta";
    std::string_view genome = std::string_view(genome_name);
    std::vector<int> sequence_lengths_vector;
    for (int i = 0; i < no_of_sequences; i++) {
        sequence_lengths_vector.push_back(sequence_lengths[i]);
    }
    std::cout << "GENOME IS: " << genome << std::endl;
    //std::vector<std::array<int, 2>> foundInGenome
    //        = graph->read_mapping_pantools(read, genome, sequence_lengths_vector);
    std::vector<std::array<int, 2>> foundInGenome = graph->read_mapping_pantools_both_sides(read, genome, sequence_lengths_vector);
    for (int i = 0; i < (int)foundInGenome.size(); i++) {
    }
    std::cout << "Readmapping DONE!" << std::endl;
    //results = new HitsPerSequence[no_of_sequences];
    std::vector<int> tmp_hits;
    std::vector<int> seq_offsets;
    //*results = new HitsPerSequence();
    int* offsets = new int[no_of_sequences];
    for (int i = 0; i < no_of_sequences; i++) {
        std::cout << sequence_lengths[i] << std::endl;
        int counter = 0;
        for (int j = 0; j < (int)foundInGenome.size(); j++) {
            std::cout << "SEQUENCE: " << foundInGenome[j][0] << std::endl;
            if (foundInGenome[j][0] == i + 1) {
                tmp_hits.push_back(foundInGenome[j][1]);
                counter++;
                std::cout << "Sequence: " << foundInGenome[j][0]
                          << ", Position: " << foundInGenome[j][1] << std::endl;
            }
        }
        seq_offsets.push_back(counter);
    }
        int* arraymarray = new int[tmp_hits.size()];
        std::copy(tmp_hits.begin(), tmp_hits.end(), arraymarray);
        std::copy(seq_offsets.begin(), seq_offsets.end(), offsets);
        results[0].hitsPerSequence = arraymarray;
        results[0].sequenceOffsets = offsets;

        //'const int no_of_hits = foundInGenome.size();
        //const int const_size = 2;
        /*int hits_in_seq[const_size][no_of_hits] = {0};
        for (int j = 0; j < (int)tmp_hits.size(); j++) {
            hits_in_seq[j][1] = foundInGenome[j][0];
        }*/

}

std::shared_ptr<AnnotatedDBG> load_dbg() {

    int argc = 9;
    char** argv = (char**)malloc(argc * sizeof(const char*));
    argv[0] = (char *) "/Users/patrick_flege/git/metagraph/metagraph/cmake-build-debug/metagraph_DNA5";
    argv[1] = (char *)"query";
    argv[2] = (char *)"--query-mode";
    argv[3] = (char *)"coords";
    argv[4] = (char *)"-i";

    /*
    argv[5] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/chloroplasts_changed_data_17_09_2025/graph.dbg";
    argv[6] = (char *)"-a";
    argv[7] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/chloroplasts_changed_data_17_09_2025/anno.column_coord.annodbg";
    argv[8] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/chloroplasts_changed_data_17_09_2025/test.fasta";
    auto config = std::make_unique<mtg::cli::Config>(argc, argv);
    std::string filename = "/Users/patrick_flege/git/patrick-pan-tools/chloroplasts_changed_data_17_09_2025/graph.dbg";

    argv[5] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/succinct_data/graph.dbg";
    argv[6] = (char *)"-a";
    argv[7] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/succinct_data/anno.column_coord.annodbg";
    argv[8] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/succinct_data/test.fasta";
    auto config = std::make_unique<mtg::cli::Config>(argc, argv);
    std::string filename = "/Users/patrick_flege/git/patrick-pan-tools/succinct_data/graph.dbg";
    */

    argv[5] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/pecto_test_dir/graph.dbg";
    argv[6] = (char *)"-a";
    argv[7] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/pecto_test_dir/anno.brwt_coord.annodbg";
    argv[8] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/pecto_test_dir/test.fasta";
    auto config = std::make_unique<mtg::cli::Config>(argc, argv);
    std::string filename = "/Users/patrick_flege/git/patrick-pan-tools/pecto_test_dir/graph.dbg";

    std::shared_ptr<DBGSuccinct> boss_graph = mtg::cli::load_critical_graph_from_file<DBGSuccinct>(config->infbase);
    std::shared_ptr<DeBruijnGraph> dbg = mtg::cli::load_critical_dbg(filename);
    std::shared_ptr<AnnotatedDBG> anno_graph =  mtg::cli::initialize_annotated_dbg(dbg, *config);

    std::cout << anno_graph->label_exists("pecto_dickeya_input/genomes/GCF_000803215.1_ASM80321v1_genomic.fna") << std::endl;
    std::cout << "Loading annotated graph from " << filename << std::endl;
    uint numba = anno_graph -> get_graph().num_nodes();
        //dbg.get->num_nodes();
    std::cout << "Number of nodes: " << numba << std::endl;
    //DBGSuccinct boss = *boss_graph;
    //return DBGSuccinct(boss_graph.get(), DeBruijnGraph::BASIC);
    //anno_graph->get_sequence_for_coords("FOOO", 0,2);
    return anno_graph;
}