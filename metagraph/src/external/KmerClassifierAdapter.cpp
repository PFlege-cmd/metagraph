//
// Created by Patrick Flege on 22/04/2026.
//

#include "KmerClassifierAdapter.h"

#include "KmerClassifier.hpp"
#include "cli/DeBruijnGraphWrapper.h"
#include "cli/graph_glue.hpp"

// struct TotalKmerMatrix {
//     int* core;
//     int* accessory;
//     int* unique;
//     explicit TotalKmerMatrix(int* core_kmers, int* accessory_kmers, int* unique_kmers) :
//     core(core_kmers), accessory(accessory_kmers), unique(unique_kmers) {
//     };
// };
//
// struct DistinctKmerMatrix {
//     int* core;
//     int* accessory;
//     int* unique;
//     explicit DistinctKmerMatrix(int* core_kmers, int* accessory_kmers, int* unique_kmers) :
//     core(core_kmers), accessory(accessory_kmers), unique(unique_kmers){
//     };
// };
//
// struct KmerMatrix {
//     DistinctKmerMatrix distinctKmerMatrix;
//     TotalKmerMatrix totalKmerMatrix;
//
//     explicit KmerMatrix(const DistinctKmerMatrix* distinctKmerMatrix_ptr,
//                         const TotalKmerMatrix* total_kmer_matrix_ptr)
//         : distinctKmerMatrix(*distinctKmerMatrix_ptr), totalKmerMatrix(*total_kmer_matrix_ptr) {
//     };
// };

struct KmerMatrix {
    int* core_total;
    int* accessory_total;
    int* unique_total;

    int* core_distinct;
    int* accessory_distinct;
    int* unique_distinct;

    int* distinct_shared;
    int* distinct_total;
    int* all_shared;
    int* all_total;

    // this->distinct_shared = std::vector<std::vector<int>>(this->genome_number, std::vector<int>(this->genome_number, 0));
    // this->distinct_total = std::vector<std::vector<int>>(this->genome_number, std::vector<int>(this->genome_number, 0));
    // this->all_shared_matrix = std::vector<std::vector<int>>(this->genome_number, std::vector<int>(this->genome_number, 0));
    // this->all_total_matrix = std::vector<std::vector<int>>(this->genome_number, std::vector<int>(this->genome_number, 0));


        explicit KmerMatrix(int* core_total_ptr,
                            int* accessory_total_ptr,
                            int* unique_total_ptr,
                            int* core_distinct_ptr,
                            int* accessory_distinct_ptr,
                            int* unique_distinct_ptr,
                            int* distinct_shared_ptr,
                            int* distinct_total_ptr,
                            int* all_shared_ptr,
                            int* all_total_ptr)
            : core_total(core_total_ptr),
    accessory_total(accessory_total_ptr),
    unique_total(unique_total_ptr),
    core_distinct(core_distinct_ptr),
    accessory_distinct(accessory_distinct_ptr),
    unique_distinct(unique_distinct_ptr),
    distinct_shared(distinct_shared_ptr),
    distinct_total(distinct_total_ptr),
    all_shared(all_shared_ptr),
    all_total(all_total_ptr)
    {
        };
    };

// int* distinct_shared;
// int* distinct_total;
// int* all_shared;
// int* all_total;
extern "C" {
    KmerMatrix* classify_kmers_genome(char * database_path, int genome_number) {
        std::string data_path = std::string(database_path);
        graph_glue glue = graph_glue(0, NULL);

        static std::shared_ptr<AnnotatedDBG> graph_coord = glue.load_coord_dbg(data_path);

        static std::shared_ptr<AnnotatedDBG> graph = glue.load_dbg(data_path);

        auto wrapper = DeBruijnGraphWrapper(*graph);
        auto kmerClassifier = KmerClassifier(wrapper, genome_number, 1);
        //int num_genomes = kmerClassifier.get_num_genomes();
        kmerClassifier.set_num_genomes(genome_number);
        kmerClassifier.create_kmer_classification_matrix();
        auto total_kmer_matrix = kmerClassifier.get_total_kmer_matrix();
        auto distinct_kmer_matrix = kmerClassifier.get_distinct_kmer_matrix();

        for (size_t i = 0; i < total_kmer_matrix.size(); i++) {
            for (size_t j = 0; j < total_kmer_matrix[i].size(); j++) {
                std::cout << "Kmer entry at :" << i << j << "--" << total_kmer_matrix[i][j] << std::endl;
            }
        }
        int* core_total_array = new int[total_kmer_matrix[0].size()];
        std::move(total_kmer_matrix[0].begin(), total_kmer_matrix[0].end(), core_total_array);

        int* accessory_total_array = new int[total_kmer_matrix[1].size()];
        std::move(total_kmer_matrix[1].begin(), total_kmer_matrix[1].end(), accessory_total_array);

        int* unique_total_array = new int[total_kmer_matrix[2].size()];
        std::move(total_kmer_matrix[2].begin(), total_kmer_matrix[2].end(), unique_total_array);

        int* core_distinct_array = new int[distinct_kmer_matrix[0].size()];
        std::move(distinct_kmer_matrix[0].begin(), distinct_kmer_matrix[0].end(), core_distinct_array);

        int* accessory_distinct_array = new int[distinct_kmer_matrix[1].size()];
        std::move(distinct_kmer_matrix[1].begin(), distinct_kmer_matrix[1].end(), accessory_distinct_array);

        int* unique_distinct_array = new int[distinct_kmer_matrix[2].size()];
        std::move(distinct_kmer_matrix[2].begin(), distinct_kmer_matrix[2].end(), unique_distinct_array);

        std::vector<std::vector<int>> distinct_shared = kmerClassifier.get_distinct_shared_matrix();
        std::vector<std::vector<int>> distinct_total = kmerClassifier.get_distinct_total_matrix();
        std::vector<std::vector<int>> all_shared = kmerClassifier.get_all_shared_matrix();
        std::vector<std::vector<int>> all_total = kmerClassifier.get_all_total_matrix();

        KmerMatrix *kmerMatrix = new KmerMatrix(core_total_array,
            accessory_total_array,
            unique_total_array,
            core_distinct_array,
            accessory_distinct_array,
            unique_distinct_array,
            kmerClassifier.flatten(distinct_shared),
            kmerClassifier.flatten(distinct_total),
            kmerClassifier.flatten(all_shared),
            kmerClassifier.flatten(
                        all_total)
            ); // Ahh... this new keyword creates a pointer, therefore it did not work
        return kmerMatrix;
}

    // int* distinct_shared;
    // int* distinct_total;
    // int* all_shared;
    // int* all_total;

    void free_matrix(KmerMatrix* matrix) {
        delete matrix->core_total;
        delete matrix->accessory_total;
        delete matrix->unique_total;

        delete matrix->core_distinct;
        delete matrix->accessory_distinct;
        delete matrix->unique_distinct;
    }
}