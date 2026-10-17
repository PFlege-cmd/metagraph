//
// Created by Patrick Flege on 22/04/2026.
//

#include "KmerClassifierAdapter.h"

#include "KmerClassifier.hpp"
#include "cli/DeBruijnGraphWrapper.h"
#include "cli/graph_glue.hpp"

extern "C" {
NativeStructures::KmerMatrix* classify_kmers_genome_phenotype(char * database_path, int genome_number, char** phenotypes) {
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


        NativeStructures::KmerMatrix *kmerMatrix = new NativeStructures::KmerMatrix(core_total_array,
            accessory_total_array,
            unique_total_array,
            core_distinct_array,
            accessory_distinct_array,
            unique_distinct_array,
            kmerClassifier.flatten(distinct_shared),
            kmerClassifier.flatten(distinct_total),
            kmerClassifier.flatten(all_shared),
            kmerClassifier.flatten(
                        all_total), nullptr,
                        shared_phenotype_array, specific_phenotype_array,exclusive_phenotype_array, 0, 0
            ); // Ahh... this new keyword creates a pointer, therefore it did not work
        return kmerMatrix;
    }
}
extern "C" {

NativeStructures::KmerMatrix* classify_kmers_genome(char * database_path, int genome_number) {

        std::string data_path = std::string(database_path);
        graph_glue glue = graph_glue(0, NULL);

        static std::shared_ptr<AnnotatedDBG> graph = glue.load_dbg(data_path);

        auto wrapper = DeBruijnGraphWrapper(*graph);
        auto kmerClassifier = KmerClassifier(wrapper, genome_number, 1);
        kmerClassifier.set_num_genomes(genome_number);

        size_t initial_node_idx{kmerClassifier.get_first_real_node_idx()};
        std::cout << "Initial node index: " << initial_node_idx << std::endl;
        kmerClassifier.create_kmer_classification_matrix_parallel(initial_node_idx);

        return KmerClassifierAdapter::export_matrix(kmerClassifier);
}

    void free_matrix(NativeStructures::KmerMatrix* matrix) {
        delete matrix->core_total;
        delete matrix->accessory_total;
        delete matrix->unique_total;

        delete matrix->core_distinct;
        delete matrix->accessory_distinct;
        delete matrix->unique_distinct;
    }
}