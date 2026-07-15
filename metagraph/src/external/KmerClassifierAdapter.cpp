//
// Created by Patrick Flege on 22/04/2026.
//

#include "KmerClassifierAdapter.h"

#include "KmerClassifier.hpp"
#include "cli/DeBruijnGraphWrapper.h"
#include "cli/graph_glue.hpp"
#include <fstream>


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

    int* shared_phenotypes;
    int* specific_phenotypes;
    int* exclusive_phenotypes;
        explicit KmerMatrix(int* core_total_ptr,
                            int* accessory_total_ptr,
                            int* unique_total_ptr,
                            int* core_distinct_ptr,
                            int* accessory_distinct_ptr,
                            int* unique_distinct_ptr,
                            int* distinct_shared_ptr,
                            int* distinct_total_ptr,
                            int* all_shared_ptr,
                            int* all_total_ptr,
                            int* shared_phenotypes_ptr,
                            int* specific_phenotypes_ptr,
                            int* exclusive_phenotypes_ptr)
            :
    core_total(core_total_ptr),
    accessory_total(accessory_total_ptr),
    unique_total(unique_total_ptr),
    core_distinct(core_distinct_ptr),
    accessory_distinct(accessory_distinct_ptr),
    unique_distinct(unique_distinct_ptr),
    distinct_shared(distinct_shared_ptr),
    distinct_total(distinct_total_ptr),
    all_shared(all_shared_ptr),
    all_total(all_total_ptr),
    shared_phenotypes(shared_phenotypes_ptr),
    specific_phenotypes(specific_phenotypes_ptr),
    exclusive_phenotypes(exclusive_phenotypes_ptr)
    {
        };
    };

// int* distinct_shared;
// int* distinct_total;
// int* all_shared;
// int* all_total;
extern "C" {
    KmerMatrix* classify_kmers_genome_phenotype(char * database_path, int genome_number, char** phenotypes) {
        map<std::string, std::vector<int>> phenotype_genome_map;
        map<int, std::string> phenotype_map;
        std::vector<std::string> phenotype_vector = vector<std::string>(genome_number);
        for (int i = 0; i < genome_number; i++) {
            phenotype_vector[i] = phenotypes[i];
            phenotype_map[i] = phenotypes[i];
            if (phenotype_genome_map.find(phenotypes[i]) != phenotype_genome_map.end()) {
                std::cout << "Added phenotype: " << phenotypes[i] << std::endl;
                std::cout << "At: " << i << std::endl;
                phenotype_genome_map[phenotypes[i]].push_back(i);
            } else {
                std::cout << "New vector: " << i << std::endl;
                std::cout << "Genome: " << phenotypes[i] << std::endl;
                phenotype_genome_map[phenotypes[i]] = std::vector<int>({i});
            }

        }

        map<int, std::vector<int>> num_phenotype_genomes;
        map<int, int> genome_num_phenotype_map;
        int counter = 0;
        for (auto const& [key, value] : phenotype_genome_map) {
            std::cout << key << " "  << "\n";
            std::cout << "Enumerating:" << std::endl;
            for (auto const& genome : value) {
                std::cout << genome<< " ";
                genome_num_phenotype_map[genome] = counter;
            }
            num_phenotype_genomes[counter] = std::vector<int>(value);
            std::cout << "Mapped successfully!" << "\n";
            std::cout << database_path << "\n";
            counter++;
        }

        map<int, int> phenotype_frequency_map;
        std::vector<int> phenotype_frequencies;

        for (auto const& [phenotype, genomes] : num_phenotype_genomes) {
            std::cout << "Genome: " <<  phenotype << " "  << "\n";
            std::cout << "Phenotypes:" << std::endl;
            phenotype_frequency_map[phenotype] = (int) genomes.size();
            phenotype_frequencies.push_back((int) genomes.size());
            for (auto const& phenotype: genomes) {
                std::cout << phenotype<< " " << std::endl;
            }
        }

        for (auto const& [key, value] : phenotype_map) {
            std::cout << key << " " << "\n";
            std::cout << "Phenotype:" << std::endl;
            std::cout << value << "\n";
        }

        std::string data_path = std::string(database_path);
        graph_glue glue = graph_glue(0, NULL);

        static std::shared_ptr<AnnotatedDBG> graph = glue.load_dbg(data_path);

        auto wrapper = DeBruijnGraphWrapper(*graph);
        auto kmerClassifier = KmerClassifier(wrapper, genome_number, 1);
        kmerClassifier.set_num_genomes(genome_number);
        auto matrix_pair = kmerClassifier.create_kmer_and_phenotype_classification_matrix(phenotype_frequencies, genome_num_phenotype_map); //TODO: Overload this.
        auto total_kmer_matrix = matrix_pair.first;
        auto phenotype_matrix = matrix_pair.second;
        auto distinct_kmer_matrix = kmerClassifier.get_distinct_kmer_matrix();

        //genome_num_phenotype_map
        //phenotype_frequency_map
        //phenptype_frequencies

        std::vector<std::vector<int>> distinct_shared = kmerClassifier.get_distinct_shared_matrix();
        std::vector<std::vector<int>> distinct_total = kmerClassifier.get_distinct_total_matrix();
        std::vector<std::vector<int>> all_shared = kmerClassifier.get_all_shared_matrix();
        std::vector<std::vector<int>> all_total = kmerClassifier.get_all_total_matrix();

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

        int* shared_phenotype_array = new int[phenotype_matrix[0].size()]; // Does not actually matter since all have the same size.
        std::move(phenotype_matrix[0].begin(), phenotype_matrix[0].end(), shared_phenotype_array);

        int* specific_phenotype_array = new int[phenotype_matrix[1].size()];
        std::move(phenotype_matrix[1].begin(), phenotype_matrix[1].end(), specific_phenotype_array);

        int* exclusive_phenotype_array = new int[phenotype_matrix[2].size()];
        std::move(phenotype_matrix[2].begin(), phenotype_matrix[2].end(), exclusive_phenotype_array);


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
                        all_total),
                        shared_phenotype_array, specific_phenotype_array,exclusive_phenotype_array
            ); // Ahh... this new keyword creates a pointer, therefore it did not work
        return kmerMatrix;
    }
}
extern "C" {

    KmerMatrix* classify_kmers_genome(char * database_path, int genome_number) {
        //map<std::string, std::vector<int>> kmer_map;
        //
        // std::string input_kmer;
        // //map<std::string, int> kmer_map;
        // ifstream ifs("/Users/patrick_flege/git/patrick-pan-tools/covid/bin/out-kmers.txt");
        //
        // while (getline(ifs, input_kmer)) {
        //     std::cout << input_kmer << std::endl;
        //     std::string kmer_component = input_kmer.substr(0, 13);
        //     //std::cout << "Extracted compoennet." << kmer_component << ::endl;
        //     //std::cout << "Position: " << input_kmer.find(' ') << std::endl;
        //     std::string numba = input_kmer.substr(13, input_kmer.size() - 13);
        //     //std::cout << "NUMBA: " << numba << std::endl;
        //     int count_kmer  = std::stoi(numba);
        //     //std::cout << "Extracted count." << std::endl;
        //     kmer_map[input_kmer] = count_kmer;
        //     std::cout << kmer_component << std::endl;
        //     std::cout << kmer_map[input_kmer] << std::endl;
        //     kmer_map.insert(std::pair<std::string, int>(input_kmer, count_kmer));
        // }

        // int test = kmer_map["AAAAAAAAAAAAA"];
        // if (test == 0) {
        //     std::cout << "NULLNULLNULL" << std::endl;
        //     throw runtime_error("KmerClassifier doesn't exist");
        // }
        std::string data_path = std::string(database_path);
        graph_glue glue = graph_glue(0, NULL);

        //static std::shared_ptr<AnnotatedDBG> graph_coord = glue.load_coord_dbg(data_path);

        static std::shared_ptr<AnnotatedDBG> graph = glue.load_dbg(data_path);

        auto wrapper = DeBruijnGraphWrapper(*graph);
        auto kmerClassifier = KmerClassifier(wrapper, genome_number, 1);
        //int num_genomes = kmerClassifier.get_num_genomes();
        kmerClassifier.set_num_genomes(genome_number);
        //kmerClassifier.set_kmer_map(kmer_map);
        kmerClassifier.create_kmer_classification_matrix(); //TODO: Overload this.
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
                        all_total), nullptr, nullptr, nullptr
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