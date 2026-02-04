//
// Created by Patrick Flege on 23/07/2025.
//

#include "load/load_annotated_graph.hpp"
#include <filesystem>
#include <iostream>

#include "load/load_graph.hpp"
#include "graph_glue.hpp"

#include "DeBruijnGraphWrapper.h"
#include "SequenceRetriever.h"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

#include <spdlog/sinks/stdout_color_sinks.h>

#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "cli/config/config.hpp"
#include "cli/build.hpp"
#include "cli/annotate.hpp"
#include "cli/stats.hpp"
#include "cli/augment.hpp"
#include "cli/clean.hpp"
#include "cli/merge.hpp"
#include "cli/align.hpp"
#include "cli/query.hpp"
#include "cli/assemble.hpp"
#include "cli/server.hpp"
#include "cli/transform_graph.hpp"
#include "cli/transform_annotation.hpp"
#include <chrono>
#include <ranges>
using namespace mtg::graph;
namespace fs = std::filesystem;
using namespace std;
using node_index = SequenceGraph::node_index;
extern "C"{
    unsigned long get_maximum_kmer_frequency(char * database_path) {
        std::cout << "get_maximum_kmer_frequency in C++!" << std::endl;
        std::string data_path = std::string(database_path);
        graph_glue glue = graph_glue(0, NULL);

        unsigned long long num_top_labels = 4294967295;
        const double discovery_fraction = 0.699999999999996;
        const double presence_fraction = 0.0;

        static std::shared_ptr<AnnotatedDBG> graph = glue.load_dbg(data_path);
        uint64 n = graph->get_graph().num_nodes();
        uint64 current = 1;
        long max_freq = 0;
        std::cout <<"Reaches line 53!" << std::endl;

        while (current < n) {
            long sum_total = 0;
            auto single_node_vec = std::vector<node_index>{current};
            auto count_vector = graph->get_kmer_counts(single_node_vec, num_top_labels, discovery_fraction,
                                                  presence_fraction);
            for (int i = 0; i < (int) count_vector.size(); i++) {
                auto stuff = std::get<2>(count_vector[i]);
                for (int j = 0; j < (int) stuff.size(); j++) {
                    sum_total += stuff[j];
                }
            }
            max_freq = sum_total > max_freq ? sum_total : max_freq;
            current += 1;
        }
        return max_freq;
    }
}

extern "C"{
    __attribute__((visibility("default")))
    __attribute__((used))
    int get_very_special_secret_message() {
        graph_glue glue = graph_glue(0, NULL);
        glue.load_dbg();
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
        graph_glue glue = graph_glue(0, NULL);
        static std::shared_ptr<AnnotatedDBG> graph = glue.load_dbg();
        HitsPerSequence* results = new HitsPerSequence();
        glue.do_pantools_work(genome_name, sequence_lengths, no_of_sequences, results, read, graph);

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

    const char * retrieve_sequence_for_coordinates_with_anchor(char * anchor_sequence, long anchor_position, char * genome_name, long start, long end) {
        graph_glue glue = graph_glue(0, NULL);
        static std::shared_ptr<AnnotatedDBG> graph = glue.load_dbg();
        AnnotatedDBG* graph_ptr = graph.get();
        DeBruijnGraphWrapper wrapper(*graph_ptr);// =  new DeBruijnGraphWrapper(*graph_ptr);
        SequenceRetriever retriever(wrapper);
        retriever.setGenome(genome_name);
        retriever.setAnchorPosition(anchor_position);
        retriever.setSequenceStart(start);
        retriever.setSequenceEnd(end);

        uint64_t start_id = retriever.retrieveAnchorId(anchor_sequence);
        uint64_t end_id = retriever.retrieveStartId(start_id);

        //std::cout << "END_ID: " << end_id << std::endl;
        std::string seq_end_id = graph_ptr->get_graph().get_node_sequence(end_id);
        //std::cout << "Sequence end ID: " << seq_end_id << std::endl;

        std::string retrieved_sequence = retriever.retrieveRegionForRange(end_id);
        std::cout << "Retrieved sequence: " << retrieved_sequence << std::endl;

        char * c_str_retrieved_sequence = (char*)malloc(retrieved_sequence.length() + 1);
        std::strcpy(c_str_retrieved_sequence, retrieved_sequence.c_str());
        return c_str_retrieved_sequence;
    }

}





extern "C"{
    __attribute__((visibility("default")))
    __attribute__((used))
    const char *  retrieve_sequence_for_coordinates(char * genome_name, long start, long end) {
        graph_glue glue = graph_glue(0, NULL);
        static std::shared_ptr<AnnotatedDBG> graph = glue.load_dbg();
        std::string seq = graph->get_sequence_for_coords(genome_name, start, end);
        std::cout << "Checking genome: " << std::endl;
        std::cout << genome_name << std::endl;
        char * result = (char*)malloc(seq.length() + 1);
        std::strcpy(result, seq.c_str());
        return result;
    }
}

extern "C"{
    void construct_graph_from_pantools() {

        /*

        build
-v
-p
4
-k
19
-o
graph
../pecto_dickeya_input/genomes/GCF_000147055.1_ASM14705v1_genomic.fasta
../pecto_dickeya_input/genomes/GCF_000803215.1_ASM80321v1_genomic.fasta
../pecto_dickeya_input/genomes/GCF_000808115.1_ASM80811v1_genomic.fasta
../pecto_dickeya_input/genomes/GCF_000808375.1_ASM80837v1_genomic.fasta
../pecto_dickeya_input/genomes/GCF_002904195.1_ASM290419v1_genomic.fasta
../pecto_dickeya_input/genomes/GCF_003595035.1_ASM359503v1_genomic.fasta
         */
        int argc = 14;
        char** argv = (char**)malloc(argc * sizeof(const char*));
        argv[0] = (char *) "/Users/patrick_flege/git/metagraph/metagraph/cmake-build-debug/metagraph_DNA5";
        argv[1] = (char *)"build";
        argv[2] = (char *)"-v";
        argv[3] = (char *)"-p";
        argv[4] = (char *)"4";
        argv[5] = (char *) "-k";
        argv[6] = (char *)"19";
        argv[7] = (char *)"-o";
        argv[8] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_output/graph_ara";
        argv[9] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_input/GCA_028009825.2_Col-CC_genomic.fna";
        argv[10] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_input/GCA_051624255.1_T8_assembly_genomic.fna";
        argv[11] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_input/GCA_051624265.1_F8_assembly_genomic.fna";
        argv[12] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_input/GCA_946409825.1_Tanz-1.10024.PacbioHiFiAssembly_genomic.fna";
        argv[13] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_input/GCF_000001735.4_TAIR10.1_genomic.fna";


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

        //argv[5] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/pecto_test_dir/graph.dbg";
        //argv[6] = (char *)"-a";
        //argv[7] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/pecto_test_dir/anno.brwt_coord.annodbg";
        //argv[8] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/pecto_test_dir/test.fasta";
        auto config = std::make_unique<mtg::cli::Config>(argc, argv);
        auto beg = std::chrono::high_resolution_clock::now();
        build_graph(config.get());
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - beg);
        std::cout << "Graph CONSTRUCTION TIME: " << duration.count() << " microseconds" << std::endl;
    }
}

extern "C"{
    void create_graph_external(int cmd_arg_number, char** cmd_args) {
        graph_glue glue = graph_glue(cmd_arg_number, cmd_args);

        CommandLineInterfaceCaller caller = CommandLineInterfaceCaller();
        AbstractCommandLineInterface * caller_interface = &caller;
        glue.set_cli_caller(*caller_interface);


        char* arg_app = (char*) fs::current_path().parent_path().append("metagraph").append("metagraph").append("cmake-build-debug").append("metagraph_DNA5").c_str();
        std::cout <<"Testing getting path metagraph:!" << std::endl;
        std::cout << arg_app << std::endl;


        std::unique_ptr<Config> config = glue.create_config();
        std::cout << "CALL SUCCESSFUL! " << std::endl;
        auto beg = std::chrono::high_resolution_clock::now();
        glue.call_cmdline_flow(config);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - beg);
        std::cout << "Graph ANNOTATION TIME FROM PANTOOLS: " << duration.count() << " microseconds" << std::endl;
    }
}

extern "C"{
    __attribute__((visibility("default")))
    void annotate_from_pantools() {
        int argc = 14;

        char** argv = (char**)malloc(argc * sizeof(const char*));
        argv[0] = (char *)"/Users/patrick_flege/git/metagraph/metagraph/cmake-build-debug/metagraph_DNA5";
        argv[1] = (char *)"annotate";
        argv[2] = (char *)"-v";
        argv[3] = (char *)"-i";
        argv[4] = (char *) "/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_output/graph_ara.dbg";
        argv[5] = (char *)"--anno-filename";
        argv[6] = (char *)"--coordinates";
        argv[7] = (char *)"-o";
        argv[8] = (char *) "/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_output/annotation";
        argv[9] = (char *) "/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_input/GCA_028009825.2_Col-CC_genomic.fna";
        argv[10] = (char *) "/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_input/GCA_051624255.1_T8_assembly_genomic.fna";
        argv[11] = (char *) "/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_input/GCA_051624265.1_F8_assembly_genomic.fna";
        argv[12] = (char *) "/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_input/GCA_946409825.1_Tanz-1.10024.PacbioHiFiAssembly_genomic.fna";
        argv[13] = (char *) "/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_input/GCF_000001735.4_TAIR10.1_genomic.fna";
        auto config = std::make_unique<mtg::cli::Config>(argc, argv);
        auto beg = std::chrono::high_resolution_clock::now();
        annotate_graph(config.get());
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - beg);
        std::cout << "Graph ANNOTATION TIME FROM PANTOOLS: " << duration.count() << " microseconds" << std::endl;
    }
}

extern "C"{
    __attribute__((visibility("default")))
        void compress_annotation_with_pantools() {

        int argc = 12;
        char** argv = (char**)malloc(argc * sizeof(const char*));

        argv[0] = (char*) "/Users/patrick_flege/git/metagraph/metagraph/cmake-build-debug/metagraph_DNA5";
        argv[1] = (char*) "transform_anno";
        argv[2] = (char*) "-v";
        argv[3] = (char*) "-p";
        argv[4] = (char*) "18";
        argv[5] = (char*) "--anno-type";
        argv[6] = (char*) "brwt_coord";
        argv[7] = (char*) "--greedy";
        argv[8] = (char*) "-o";
        argv[9] = (char*) "/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_output/anno";
        argv[10] = (char*) "--coordinates";
        argv[11] = (char*) "/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_output/annotation.column.annodbg";


        auto config = std::make_unique<mtg::cli::Config>(argc, argv);
        auto beg = std::chrono::high_resolution_clock::now();
        transform_annotation(config.get());
        //std::vector<string> fname_vec = {"/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_output/annotation.column.annodbg"};
        //config->fnames = fname_vec;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - beg);
        std::cout << "Annotation COMPRESSION TIME: " << duration.count() << " microseconds" << std::endl;
    }
}

extern "C"{
    const char *  retrieve_kmer_for_coordinates(int genome, int start) {
        std::cout << "Checking genome: " << std::endl;
        std::cout << genome << std::endl;
        std::cout << start << std::endl;
        graph_glue glue(0, NULL);
        static std::shared_ptr<AnnotatedDBG> graph = glue.load_dbg();
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

graph_glue::graph_glue(int argcount, char** argv) {
    //std::cout << "Creating Graph glue!" << std::endl;
    cmd_arguments = (char**)calloc(argcount, sizeof(const char*));
    argc = argcount;
    for (int i = 0; i < argc; i++) {
        cmd_arguments[i] = strdup(argv[i]);
        std::cout << "Arg " << i << " is: " << cmd_arguments[i] << std::endl;
    }
    //int external_arg_counter = 0;
    //for (int i = 0; argv[i] != nullptr; i++) {
    //    external_arg_counter++;
    //}
    //assert(external_arg_counter == argcount);
}

graph_glue::~graph_glue() {
    std::cout << "Destroying Graph glue!" << std::endl;

    for (int i = 0; i < argc; i++) {
        free(cmd_arguments[i]);
        cmd_arguments[i] = nullptr; //to avoid dangling pointer
    }
    free(cmd_arguments);
    cmd_arguments = nullptr;
}

void graph_glue::set_cli_caller(AbstractCommandLineInterface& cli) {
    cli_caller = &cli;
}


void graph_glue::do_pantools_work(char* genome_name,
             int* sequence_lengths,
             int no_of_sequences,
             HitsPerSequence*& results,
                      const char* read_ptr, const std::shared_ptr<AnnotatedDBG>& graph) {
    //static std::shared_ptr<AnnotatedDBG> graph = load_dbg();
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

std::shared_ptr<AnnotatedDBG> graph_glue::load_dbg() {
    int argc = 9;
    char** argv = (char**) malloc(argc * sizeof(const char*));
    //char* arg_app = (char*) fs::current_path().parent_path().append("metagraph").append("metagraph").append("cmake-build-debug").append("metagraph_DNA5").c_str();
    char* arg_app = (char*) fs::current_path().parent_path().append("metagraph").append("metagraph").append("build").append("metagraph_DNA5").c_str();

    std::cout << arg_app << std::endl;
    argv[0] = arg_app;
    //argv[0] = (char *) "/Users/patrick_flege/git/metagraph/metagraph/cmake-build-debug/metagraph_DNA5";
    argv[1] = (char*)"query";
    argv[2] = (char*)"--query-mode";
    argv[3] = (char*)"coords";
    argv[4] = (char*)"-i";

    // LUSTRE, SERVER
    /*
    argv[5] = (char *)"/lustre/BIF/nobackup/flege001/patrick-pan-tools/chloroplast_DB/graph.dbg";
    argv[6] = (char*)"-a";
    argv[7] = (char *)"/lustre/BIF/nobackup/flege001/patrick-pan-tools/chloroplast_DB/anno.column_coord.annodbg";
    argv[8] = (char *)"/lustre/BIF/nobackup/flege001/patrick-pan-tools/chloroplast_DB/test.fasta";
    */

    // LOCAL, on MAC:
     argv[5] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_20_DB/graph.dbg";
     argv[6] = (char*) "-a";
     argv[7] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_20_DB/transformed.row_diff_brwt_coord.annodbg";
     argv[8] = (char *)"/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_20_DB/test.fasta";
     auto config = std::make_unique<mtg::cli::Config>(argc, argv);
     // std::string filename
     //        = "/lustre/BIF/nobackup/flege001/patrick-pan-tools/chloroplast_DB/graph.dbg";

    std::string filename = "/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_20_DB/graph.dbg";

    std::shared_ptr<DBGSuccinct> boss_graph = mtg::cli::load_critical_graph_from_file<DBGSuccinct>(config->infbase);
    std::shared_ptr<DeBruijnGraph> dbg = mtg::cli::load_critical_dbg(filename);
    std::shared_ptr<AnnotatedDBG> anno_graph = mtg::cli::initialize_annotated_dbg(dbg, *config);
//
     return anno_graph;
 }

std::shared_ptr<AnnotatedDBG> graph_glue::load_dbg(std::string database_path) {
    std::cout << "Loading database from " << database_path << std::endl;
    int argc = 9;
    char** argv = (char**)malloc(argc * sizeof(const char*));
    //char* arg_app = (char*) fs::current_path().parent_path().append("metagraph").append("metagraph").append("cmake-build-debug").append("metagraph_DNA5").c_str();
    char* arg_app = (char*) fs::current_path().parent_path().append("metagraph").append("metagraph").append("build").append("metagraph_DNA5").c_str();

    std::string filename_local = fs::current_path().append(database_path).append("graph.dbg");

    std::cout << arg_app << "Application: "<< std::endl;
    argv[0] = arg_app;
    //argv[0] = (char *) "/Users/patrick_flege/git/metagraph/metagraph/cmake-build-debug/metagraph_DNA5";
    argv[1] = (char*)"query";
    argv[2] = (char*)"--query-mode";
    argv[3] = (char*)"coords";
    argv[4] = (char*)"-i";

    // LUSTRE, SERVER
    /*
    argv[5] = (char *)"/lustre/BIF/nobackup/flege001/patrick-pan-tools/chloroplast_DB/graph.dbg";
    argv[6] = (char*)"-a";
    argv[7] = (char *)"/lustre/BIF/nobackup/flege001/patrick-pan-tools/chloroplast_DB/anno.column_coord.annodbg";
    argv[8] = (char *)"/lustre/BIF/nobackup/flege001/patrick-pan-tools/chloroplast_DB/test.fasta";
    */

    // LOCAL, on MAC:
     //auto file_name_cstr = filename_local.c_str();
     argv[5] = (char *) "/lustre/BIF/nobackup/flege001/patrick-pan-tools/chloroplast_DB/graph.dbg";
     argv[6] = (char*) "-a";
     argv[7] = (char *) fs::current_path().append(database_path).append("count_brwt.int_brwt.annodbg").c_str();
     argv[8] = (char*) "/lustre/BIF/nobackup/flege001/patrick-pan-tools/chloroplast_DB/test.fasta";
     auto config = std::make_unique<mtg::cli::Config>(argc, argv);
     // std::string filename
     //        = "/lustre/BIF/nobackup/flege001/patrick-pan-tools/chloroplast_DB/graph.dbg";

    //std::string filename = "/Users/patrick_flege/git/patrick-pan-tools/a_thaliana_20_DB/graph.dbg";
    std::cout << config->infbase << std::endl;
    std::cout << "PREVIOUS WAS CONFIG" << std::endl;

    std::shared_ptr<DBGSuccinct> boss_graph = mtg::cli::load_critical_graph_from_file<DBGSuccinct>(config->infbase);
    std::shared_ptr<DeBruijnGraph> dbg = mtg::cli::load_critical_dbg(filename_local);
    std::shared_ptr<AnnotatedDBG> anno_graph = mtg::cli::initialize_annotated_dbg(dbg, *config);
    std::cout << "AM HERE!" << std::endl;
//
     return anno_graph;
 }

 char** graph_glue::get_cmd_arguments() {
     return cmd_arguments;
 }

 void graph_glue::set_cmd_arguments(char *arguments[]) {
     for (int i = 0; arguments[i] != NULL; i++) {
         cmd_arguments[i] = arguments[i];
     }
}

int graph_glue::get_cmd_arg_count() {
    return argc;
}

unique_ptr<mtg::cli::Config> graph_glue::create_config() {
    return cli_caller->get_config(argc, cmd_arguments);
}

void graph_glue::call_cmdline_flow(const std::unique_ptr<mtg::cli::Config>& config) {
     cli_caller->call_flow(config);
 }