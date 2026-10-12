//
// Created by Patrick Flege on 16/03/2026.
//

#pragma once
#include "PangenomeStructureAdapter.h"
#include "PhenotypeClass.h"
#include "cli/GraphWrapper.hpp"
#include "external/KmerType.h"

#include <regex>
#include <utility>
#include <vector>

using node_index_kmer = GraphWrapper::node_index_kmer;
using kmer_frequencies = GraphWrapper::kmer_frequencies;
class KmerClassifier {
public:
    typedef std::vector<std::pair<int, int>> genomes_and_frequencies;
    explicit KmerClassifier(GraphWrapper& graph_ptr);
    KmerClassifier(GraphWrapper& graph_ptr, int core_threshold, int unique_threshold);
    KmerClassifier(GraphWrapper& graph_ptr, int core_threshold, int unique_threshold, int genome_number);
    unsigned long long get_node_number();
    kmer_frequencies get_kmer_frequencies(std::vector<node_index_kmer>& nodes);
    int get_num_genomes();
    void set_num_genomes(int genome_number);
    int extractGenomeNumber(std::string genome_contig_path);


    std::vector<int> count_kmer_per_genome(node_index_kmer idx) {
        std::vector<int> kmer_counts(this->genome_number, { 0 });
        auto nodes = std::vector<node_index_kmer>({ idx });
        auto frequencies = this->get_kmer_frequencies(nodes);
        // std::cerr << "idx: " << idx << "frequencies.size() " << frequencies.size() << std::endl;
        for (size_t i = 0; i < frequencies.size(); i++) {
            auto genome_name = std::get<0>(frequencies[i]);
            auto genome_no = this->extractGenomeNumber(genome_name);
            if (genome_no == -1)
                throw std::invalid_argument("Invalid genome number");
            auto genome_count = std::get<2>(frequencies[i])[0];
            int genome_count_int = (int)genome_count;
            kmer_counts[genome_no] = genome_count_int;
        }
        return kmer_counts;
    }

    size_t get_first_real_node_idx() {
        size_t idx { 1 };
        std::string kmer_string { graph_ptr.get_node_sequence(idx) };
        while (kmer_string.find('$') < kmer_string.length()) {
            idx++;
            kmer_string = graph_ptr.get_node_sequence(idx);
        }
        return idx;
    }


    /**
     *@brief Does the first run of the pangenome-structure flow, filling in also the
     * bitset with indices, updating the occurrences of kmers already in the loop-collector for the case of n genomes.
     *
     *
     * @param loop_collector Main object of flow, carries vectors of vectors, which have structure:
     * Number of loops x genome_number - 1. Three of those lists are carried, one for core, accessory, and unique kmers.
     * @param all_node_ids Vector of all node-ids.
     * @returns loop-collector object with updated bit-frequencies in a single vector.
     */
    PangenomeStructureAdapter::PangenomeStructureCollector& initialize_kmer_frequencies(PangenomeStructureAdapter::PangenomeStructureCollector& loop_collector, const std::vector<unsigned long long>& all_node_ids);



    std::vector<int> count_kmer_per_genome_all(node_index_kmer idx) {

        std::vector<int> kmer_counts(this->genome_number, {0});
        auto nodes = std::vector<node_index_kmer>(15000);
        std::iota(nodes.begin(), nodes.end(), 15);
        //auto nodes = std::vector<node_index_kmer>({idx});
        auto frequencies = this->get_kmer_frequencies(nodes);
        for (size_t i = 0; i < frequencies.size(); i++) {
            auto genome_name =  std::get<0>(frequencies[i]);
            auto genome_no = this->extractGenomeNumber(genome_name);
            if (genome_no == -1)
                throw std::invalid_argument("Invalid genome number");
            auto genome_count = std::get<2>(frequencies[i])[0];
            int genome_count_int = (int)genome_count;
            kmer_counts[genome_no] = genome_count_int;
        }
        return kmer_counts;
    };

    std::vector<int> count_kmer_per_genome_test(node_index_kmer idx) {
        std::vector<int> kmer_counts(this->genome_number, {0});
        auto nodes = std::vector<node_index_kmer>({27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63});
        auto frequencies = this->get_kmer_frequencies(nodes);
        for (size_t i = 0; i < frequencies.size(); i++) {
            auto genome_name =  std::get<0>(frequencies[i]);
            auto genome_no = this->extractGenomeNumber(genome_name);
            if (genome_no == -1)
                throw std::invalid_argument("Invalid genome number");
            auto genome_count = std::get<2>(frequencies[i])[0];
            int genome_count_int = (int)genome_count;
            kmer_counts[genome_no] = genome_count_int;
        }
        return kmer_counts;
    };

    std::tuple<std::vector<int>, std::pair<PhenotypeClass, std::vector<int>>>
    classify_phenotype(std::vector<int>& pangenome_phenotype_frequencies,
                       std::vector<int>& kmer_counts,
                       int zeros,
                       std::vector<int>& freq_equal_indices,
                       std::vector<int>& not_zero_indices) {
        if (freq_equal_indices.size() == 1
            && zeros == (int)(pangenome_phenotype_frequencies.size() - 1)) {
            const auto count_and_phenotype_tuple
                    = std::make_tuple<std::vector<int>, std::pair<PhenotypeClass, int>>(
                            std::move(kmer_counts),
                            std::make_pair(PhenotypeClass::SPECIFIC, freq_equal_indices[0]));
        }
        if (freq_equal_indices.size() > 0) { // Shared
            const auto count_and_phenotype_tuple
                    = std::make_tuple<std::vector<int>, std::pair<PhenotypeClass, std::vector<int>>>(
                            std::move(kmer_counts),
                            std::make_pair<PhenotypeClass, std::vector<int>>(
                                    PhenotypeClass::SHARED, std::move(freq_equal_indices)));
            return count_and_phenotype_tuple;
        }
        if (zeros == (int)(pangenome_phenotype_frequencies.size() - 1)) {
            // EXCLUSIVE
            const auto count_and_phenotype_tuple
                    = std::make_tuple<std::vector<int>, std::pair<PhenotypeClass, std::vector<int>>>(
                            std::move(kmer_counts),
                            std::make_pair(PhenotypeClass::EXCLUSIVE, not_zero_indices));
            return count_and_phenotype_tuple;
        }

        return std::make_tuple<std::vector<int>, std::pair<PhenotypeClass, std::vector<int>>>(
                {}, std::make_pair(PhenotypeClass::NOT_PRESENT, std::vector<int>(0)));
    }
    std::tuple<std::vector<int>, std::pair<PhenotypeClass, std::vector<int>>> count_kmer_per_genome_and_classify_phenotype(node_index_kmer idx, std::vector<int> pangenome_phenotype_frequencies, std::map<int, int> genome_to_phenotype_map) {
        std::vector<int> kmer_counts(this->genome_number, {0});
        //std::vector<int> kmer_counts(genome_number, 0);
        auto nodes = std::vector<node_index_kmer>({idx});
        auto frequencies = this->get_kmer_frequencies(nodes);
        auto indices = std::vector<int>(pangenome_phenotype_frequencies.size());
        std::iota(std::begin(indices), std::end(indices), 0);


        std::vector<int> kmer_phenotype_occurrences = std::vector<int>(pangenome_phenotype_frequencies.size(), {0});
        for (size_t i = 0; i < frequencies.size(); i++) {
            auto genome_name =  std::get<0>(frequencies[i]);
            auto genome_no = this->extractGenomeNumber(genome_name);
            if (genome_no == -1)
                throw std::invalid_argument("Invalid genome number");

            auto genome_count = std::get<2>(frequencies[i])[0];
            int genome_count_int = (int)genome_count;
            if (genome_count_int >  0 ) {
                kmer_phenotype_occurrences[genome_to_phenotype_map[genome_no]]++;
            }
            kmer_counts[genome_no] = genome_count_int;
        }

        int zeros = 0;
        std::vector<int> freq_equal_indices = std::vector<int>();
        std::vector<int> not_zero_indices = std::vector<int>();
        std::for_each(indices.begin(),
            indices.end(),
            [
                kmer_pheno = kmer_phenotype_occurrences,
                pan_pheno = pangenome_phenotype_frequencies,
                z = &zeros,
                &freq_equal_indices,
                &not_zero_indices
                ](const auto& idx) {
                if (kmer_pheno[idx] == 0) {
                    (*z)++;
                } else {
                    if (kmer_pheno[idx] == pan_pheno[idx]) {
                        freq_equal_indices.push_back(idx);
                    } else {
                        not_zero_indices.push_back(idx);
                    }
                }
            });
        return classify_phenotype(pangenome_phenotype_frequencies,
            kmer_counts,
            zeros,
            freq_equal_indices,
            not_zero_indices); // Default: case phenotype exclusive
    }
    /**
     * @brief Classifies kmers into core, accessory or unique, depending on their occurrences.
     * @param kmer_counts vector containing kmers-frequencies for each genome.
     * @return tuple of kmertype, and pair of genome_idx and associated frequencies.
     */
    std::tuple<KmerType, genomes_and_frequencies> classify_kmers(std::vector<int> kmer_counts) {
        KmerType type;
        genomes_and_frequencies genome_frequencies({});
        int number_of_occuring_genomes = 0;
        for (size_t i = 0; i < kmer_counts.size(); i++) {
            if (kmer_counts[i] > 0){
                number_of_occuring_genomes++;
                auto occurrence_pair = std::make_pair(i, kmer_counts[i]);
                genome_frequencies.push_back(occurrence_pair);
            }
        }

        if (number_of_occuring_genomes == this->unique_threshold) {
            //std::cout << "Working" << std::endl;
            type = KmerType::UNIQUE;
        }
        else if (number_of_occuring_genomes == this->core_threshold) {
            //std::cout << "Not unique" << std::endl;
            type = KmerType::CORE;
        }
        else {
            type = KmerType::ACCESSORY;
        }

        return std::tuple<KmerType, genomes_and_frequencies>(type, genome_frequencies);
    };
    /**
     *
     * @param genome_frequencies genome and frequencies of a given kmer
     * @param gen_idx pair<int, int> of index,frequencies
     */
    void fill_triangular_matrices(const genomes_and_frequencies& genome_frequencies,
                                  std::pair<int, int> gen_idx);
    void fill_triangular_matrices_parallel(const genomes_and_frequencies& genome_frequencies,
                                           std::pair<int, int> outer_genome_idx,
                                           std::vector<std::vector<int>>& distinct_shared,
                                           std::vector<std::vector<int>>& all_shared_matrix);


    /**
     * @param type Kmertype (Unique, Accessory or Core)
     * @param genome_frequencies Vector containing pairs genome-number -  frequency of kmer
     */
    void fill_kmer_matrix(KmerType type,
                          genomes_and_frequencies genome_frequencies) {
        bool* is_present_in_genome = new bool[genome_number];

        switch (type) {
            case KmerType::CORE: {
                //int counter = 0;
                //std::cout << "Core" << std::endl;
                for_each(genome_frequencies.begin(), genome_frequencies.end(),
                         [this, genome_frequencies, is_present_in_genome](auto gen_idx) {
                             //std::cout << "CORE, at index: " << gen_idx.first << std::endl;
                             //std::cout << "CORE, adding: " << gen_idx.second << std::endl;
                             //counter =+  gen_idx.second;
                             this->total_kmer_matrix[0][gen_idx.first] += gen_idx.second;
                             this->distinct_kmer_matrix[0][gen_idx.first] += 1;
                             this->fill_triangular_matrices(genome_frequencies, gen_idx);
                            is_present_in_genome[gen_idx.first] = true;
                         });
                break;
            }
            case KmerType::ACCESSORY: {
                //std::cout << "Accessory" << std::endl;
                std::for_each(genome_frequencies.begin(), genome_frequencies.end(),
                              [this, genome_frequencies,is_present_in_genome](auto gen_idx) {
                                  /*std::cout << "Accessory at index: " << gen_idx.first
                                            << std::endl;
                                  std::cout << "Accessory, adding: " << gen_idx.second
                                            << std::endl;
                                            */

                                  this->total_kmer_matrix[1][gen_idx.first] += gen_idx.second;
                                  this->distinct_kmer_matrix[1][gen_idx.first] += 1;
                                  this->fill_triangular_matrices(genome_frequencies, gen_idx);
                                  is_present_in_genome[gen_idx.first] = true;

                              });
                break;
            }
            case KmerType::UNIQUE: {
                //std::cout << "Unique" << std::endl;
                //std::cout << "Unique, adding: " << genome_frequencies[0].second << std::endl;
                this->total_kmer_matrix[2][genome_frequencies[0].first] += genome_frequencies[0].second;
                this->distinct_kmer_matrix[2][genome_frequencies[0].first] += 1;

                this->all_shared_matrix[genome_frequencies[0].first][genome_frequencies[0].first] += genome_frequencies[0].second;
                this->all_total_matrix[genome_frequencies[0].first][genome_frequencies[0].first] += genome_frequencies[0].second;

                this->distinct_shared[genome_frequencies[0].first][genome_frequencies[0].first]++;
                this->distinct_total[genome_frequencies[0].first][genome_frequencies[0].first]++;
                is_present_in_genome[genome_frequencies[0].first] = true;
                break;
            }
            default:
                break;
        }
        std::vector<int> frequencies = std::vector<int>(this->genome_number);
        for_each(genome_frequencies.begin(), genome_frequencies.end(), [frequencies_ref = &frequencies](auto gen_idx) {
            (*frequencies_ref)[gen_idx.first] += gen_idx.second;
        });

        for (int i = 0; i < genome_number; i++) {
            for (int j = 0; j < genome_number; j++) {
                int highest = frequencies[i] > frequencies[j]? frequencies[i] : frequencies[j];
                if (highest == 0)
                    continue;

                ++this->distinct_total[i][j];
                this->all_total_matrix[i][j] += highest;
            }
        }
    }

    /**
     * Flattens 2D matrix of vectors to a 1D int array, which is returned.
     * @param matrix 2D vector of vector<int>
     */
    int* flatten(const std::vector<std::vector<int>>& matrix);

    /**
    Creates kmer-matrices, both including total counts and distinct counts (adds only one instead of frequency).
    Stores those in instance variables:
    total_kmer_matrix, distinct_kmer_matrix
     */
    std::vector<std::vector<int>> create_kmer_classification_matrix() {
        this->total_kmer_matrix
                = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));
        // auto total_kmer_matrix = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));
        this->distinct_kmer_matrix
                = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));

        this->distinct_shared
                = std::vector<std::vector<int>>(this->genome_number,
                                                std::vector<int>(this->genome_number, 0));
        this->distinct_total
                = std::vector<std::vector<int>>(this->genome_number,
                                                std::vector<int>(this->genome_number, 0));
        this->all_shared_matrix
                = std::vector<std::vector<int>>(this->genome_number,
                                                std::vector<int>(this->genome_number, 0));
        this->all_total_matrix
                = std::vector<std::vector<int>>(this->genome_number,
                                                std::vector<int>(this->genome_number, 0));

        auto number_of_kmers = this->get_node_number();
        std::cout << "Node number:" << number_of_kmers << std::endl;
        node_index_kmer current = 1;
        for (size_t i = 0; i < total_kmer_matrix.size(); i++) {
            for (size_t j = 0; j < total_kmer_matrix[i].size(); j++) {
                total_kmer_matrix[i][j] = 0;
            }
        }


        auto test_graph = this->graph_ptr.get_graph()->get_graph_ptr();
        //#pragma omp parallel for
        for (size_t i = 1; i <  number_of_kmers/2; i++) {

            std::cout <<"Current is:" << i << std::endl;

            auto counts_per_genomes = this->count_kmer_per_genome(current);

            auto classified_kmers = this->classify_kmers(counts_per_genomes);
            KmerType type = std::get<0>(classified_kmers);

            genomes_and_frequencies genome_frequencies = std::get<1>(classified_kmers);
            // TODO: FIll distinct counts as well.
            fill_kmer_matrix(type, genome_frequencies);

            current++;
        }

        return total_kmer_matrix;
    }

    /**
    Creates kmer-matrices, both including total counts and distinct counts (adds only one instead of frequency).
    Stores those in instance variables:
    total_kmer_matrix, distinct_kmer_matrix
     */
    std::vector<std::vector<int>> create_kmer_classification_matrix_parallel(size_t first_node_idx) {
        this->total_kmer_matrix
                = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));
        // auto total_kmer_matrix = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));
        this->distinct_kmer_matrix
                = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));

        this->distinct_shared
                = std::vector<std::vector<int>>(this->genome_number,
                                                std::vector<int>(this->genome_number, 0));
        this->distinct_total
                = std::vector<std::vector<int>>(this->genome_number,
                                                std::vector<int>(this->genome_number, 0));
        this->all_shared_matrix
                = std::vector<std::vector<int>>(this->genome_number,
                                                std::vector<int>(this->genome_number, 0));
        this->all_total_matrix
                = std::vector<std::vector<int>>(this->genome_number,
                                                std::vector<int>(this->genome_number, 0));

        auto number_of_kmers = this->get_node_number();
        std::cout << "Node number:" << number_of_kmers << std::endl;
        //node_index_kmer current = 1;
        for (size_t i = 0; i < total_kmer_matrix.size(); i++) {
            for (size_t j = 0; j < total_kmer_matrix[i].size(); j++) {
                total_kmer_matrix[i][j] = 0;
            }
        }

        auto kmer_string{graph_ptr.get_node_sequence(first_node_idx)};
        #pragma omp parallel
        {

            // init +  (number - init )/2)
            size_t limit{number_of_kmers/2 + 1};
            //size_t limit{first_node_idx + (number_of_kmers - first_node_idx)/2};

        auto local_total_kmer_matrix = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));
        auto local_distinct_kmer_matrix = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));

        auto local_distinct_shared = std::vector<std::vector<int>>(this->genome_number,std::vector<int>(this->genome_number, 0));
        auto local_distinct_total = std::vector<std::vector<int>>(this->genome_number,std::vector<int>(this->genome_number, 0));

        auto local_all_shared_matrix = std::vector<std::vector<int>>(this->genome_number,std::vector<int>(this->genome_number, 0));
        auto local_all_total_matrix = std::vector<std::vector<int>>(this->genome_number,std::vector<int>(this->genome_number, 0));

        #pragma omp for schedule(static, 3)
        for (size_t i = 1; i < limit; i++) {

            node_index_kmer current{i};
            using node_index = mtg::graph::SequenceGraph::node_index;
            std::vector<node_index> cv{current};

            auto counts_per_genomes  = this->count_kmer_per_genome(current);

            auto classified_kmers = this->classify_kmers(counts_per_genomes);
            KmerType type = std::get<0>(classified_kmers);

            genomes_and_frequencies genome_frequencies = std::get<1>(classified_kmers);
            // TODO: FIll distinct counts as well.
            fill_kmer_matrix_parallel(type, genome_frequencies,
                local_total_kmer_matrix,
                local_distinct_kmer_matrix,
                local_distinct_shared,
                local_distinct_total,
                local_all_shared_matrix,
                local_all_total_matrix);
        }
        #pragma omp critical
        {
            for (size_t r = 0; r < local_total_kmer_matrix.size(); r++) {
                for (size_t c = 0; c < local_total_kmer_matrix[r].size(); c++) {
                    this->total_kmer_matrix[r][c] += local_total_kmer_matrix[r][c];
                    this->distinct_kmer_matrix[r][c] += local_distinct_kmer_matrix[r][c];
                }
            }
            for (size_t r = 0; r < local_distinct_shared.size(); r++) {
                for (size_t c = 0; c < local_distinct_shared[r].size(); c++) {
                    this->distinct_shared[r][c] += local_distinct_shared[r][c];
                    this->distinct_total[r][c] += local_distinct_total[r][c];
                    this->all_shared_matrix[r][c] += local_all_shared_matrix[r][c];
                    this->all_total_matrix[r][c] += local_all_total_matrix[r][c];
                }
            }
        }
    }

        return total_kmer_matrix;
    }

    void fill_kmer_matrix_parallel(KmerType type, genomes_and_frequencies genome_frequencies,
        std::vector<std::vector<int>>& total_kmer_matrix,
        std::vector<std::vector<int>>& distinct_kmer_matrix,
        std::vector<std::vector<int>>& distinct_shared,
        std::vector<std::vector<int>>& distinct_total,
        std::vector<std::vector<int>>& all_shared_matrix,
        std::vector<std::vector<int>>& all_total_matrix
        ) {
        bool* is_present_in_genome = new bool[genome_number];
        switch (type) {
            case KmerType::CORE: {
                //int counter = 0;
                //std::cout << "Core" << std::endl;
                for_each(genome_frequencies.begin(), genome_frequencies.end(),
                         [&](auto gen_idx) {
                             //std::cout << "CORE, at index: " << gen_idx.first << std::endl;
                             //std::cout << "CORE, adding: " << gen_idx.second << std::endl;
                             //counter =+  gen_idx.second;
                             total_kmer_matrix[0][gen_idx.first] += gen_idx.second;
                             distinct_kmer_matrix[0][gen_idx.first] += 1;
                             this->fill_triangular_matrices_parallel(genome_frequencies, gen_idx, distinct_shared, all_shared_matrix);
                            is_present_in_genome[gen_idx.first] = true;
                         });
                break;
            }
            case KmerType::ACCESSORY: {
                //std::cout << "Accessory" << std::endl;
                std::for_each(genome_frequencies.begin(), genome_frequencies.end(),
                              [&](auto gen_idx) {
                                  // std::cout << "Accessory at index: " << gen_idx.first
                                  //           << std::endl;
                                  // std::cout << "Accessory, adding: " << gen_idx.second
                                  //           << std::endl;


                                  total_kmer_matrix[1][gen_idx.first] += gen_idx.second;
                                  distinct_kmer_matrix[1][gen_idx.first] += 1;
                                  this->fill_triangular_matrices_parallel(genome_frequencies, gen_idx, all_shared_matrix, all_total_matrix);
                                  is_present_in_genome[gen_idx.first] = true;

                              });
                break;
            }
            case KmerType::UNIQUE: {
                // std::cout << "Unique" << std::endl;
                // std::cout << "Unique, adding: " << genome_frequencies[0].second << std::endl;
                total_kmer_matrix[2][genome_frequencies[0].first] += genome_frequencies[0].second;
                distinct_kmer_matrix[2][genome_frequencies[0].first] += 1;

                all_shared_matrix[genome_frequencies[0].first][genome_frequencies[0].first] += genome_frequencies[0].second;
                all_total_matrix[genome_frequencies[0].first][genome_frequencies[0].first] += genome_frequencies[0].second;

                distinct_shared[genome_frequencies[0].first][genome_frequencies[0].first]++;
                distinct_total[genome_frequencies[0].first][genome_frequencies[0].first]++;
                is_present_in_genome[genome_frequencies[0].first] = true;
                break;
            }
            default:
                break;
        }
        std::vector<int> frequencies = std::vector<int>(this->genome_number);
        for_each(genome_frequencies.begin(), genome_frequencies.end(), [frequencies_ref = &frequencies](auto gen_idx) {
            (*frequencies_ref)[gen_idx.first] += gen_idx.second;
        });

        for (int i = 0; i < genome_number; i++) {
            for (int j = 0; j < genome_number; j++) {
                int highest = frequencies[i] > frequencies[j]? frequencies[i] : frequencies[j];
                if (highest == 0)
                    continue;

                ++distinct_total[i][j];
                all_total_matrix[i][j] += highest;
            }
        }
    }


    void fill_phenotypes_matrix(
            std::tuple_element<1, std::tuple<std::vector<int>, std::pair<PhenotypeClass, std::vector<int>>>>::type
                    & phenoclassification_and_allele) {
        PhenotypeClass pheno_type = std::get<0>(phenoclassification_and_allele);
        auto alleles = std::get<1>(phenoclassification_and_allele);
        switch (pheno_type) {
            case PhenotypeClass::SHARED:
                //std::cout << "Adding SHARED: " << std::endl;

                for (int allele : alleles) {
                    //std::cout << alleles[0] << std::endl;
                    phenotypes[0][allele]++;
                };
                break;
            case PhenotypeClass::SPECIFIC:
                std::cout << "Adding SPECIFIC: " << alleles[0] << std::endl;
                phenotypes[1][alleles[0]]++;
                break;

            case PhenotypeClass::EXCLUSIVE:
                std::cout << "Adding EXCLUSIVE: " << alleles[0] << std::endl;
                phenotypes[2][alleles[0]]++;
                break;

            default:
                std::cout << "Unknown PhenotypeClass" << std::endl;
        }
    }
    std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> create_kmer_and_phenotype_classification_matrix(std::vector<int> pangenome_phenotype_frequencies, std::map<int, int> genome_to_phenotype_map) {
        this->total_kmer_matrix = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));
        this->distinct_kmer_matrix = std::vector<std::vector<int>>(3, std::vector<int>(this->genome_number, 0));

        this->distinct_shared = std::vector<std::vector<int>>(this->genome_number, std::vector<int>(this->genome_number, 0));
        this->distinct_total = std::vector<std::vector<int>>(this->genome_number, std::vector<int>(this->genome_number, 0));
        this->all_shared_matrix = std::vector<std::vector<int>>(this->genome_number, std::vector<int>(this->genome_number, 0));
        this->all_total_matrix = std::vector<std::vector<int>>(this->genome_number, std::vector<int>(this->genome_number, 0));

        this->phenotypes = std::vector(3, std::vector<int>(pangenome_phenotype_frequencies.size(), 0));

        auto number_of_kmers = this->get_node_number();
        node_index_kmer current = 1;
        for (size_t i = 0; i < total_kmer_matrix.size(); i++) {
            for (size_t j = 0; j < total_kmer_matrix[i].size(); j++) {
                total_kmer_matrix[i][j] = 0;
            }
        }
        for (int i : pangenome_phenotype_frequencies) {
            std::cout << "Pangenome phenotype frequency: " << i << std::endl;
            std::cout << i << std::endl;
        }

        while (current < number_of_kmers) {

            auto counts_and_pheno = this->count_kmer_per_genome_and_classify_phenotype(current, pangenome_phenotype_frequencies, genome_to_phenotype_map);
            auto counts_per_genomes = std::get<0>(counts_and_pheno);
            auto phenoclassification_and_allele = std::get<1>(counts_and_pheno);

            auto classified_kmers = this->classify_kmers(counts_per_genomes);
            KmerType type = std::get<0>(classified_kmers);
            genomes_and_frequencies genome_frequencies = std::get<1>(classified_kmers);

            fill_phenotypes_matrix(phenoclassification_and_allele);
            fill_kmer_matrix(type, genome_frequencies);

            current++;
        }

        const auto matrix_pair = std::make_pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>(
                        std::move(this->total_kmer_matrix), std::move(this->phenotypes));
        return matrix_pair;
    }

    void set_kmer_map(std::map<std::string, std::vector<int>>& kmer_map) {
        std::string input_kmer;
        std::ifstream ifs("/Users/patrick_flege/git/patrick-pan-tools/covid/bin/out-kmers.txt");

        while (getline(ifs, input_kmer)) {
            std::cout << input_kmer << std::endl;
            std::string kmer_component = input_kmer.substr(0, 13);
            std::string numba = input_kmer.substr(13, input_kmer.size() - 13);
            auto vec = std::vector{std::stoi(numba)};
            this->kmer_map.insert(std::pair<std::string, std::vector<int>>(kmer_component, vec));

            std::cout << kmer_component << std::endl;
            std::cout <<  this->kmer_map[kmer_component][0] << std::endl;
        }
    }

    std::map<std::string, std::vector<int>>& get_kmer_map() {
        return this->kmer_map;
    }

    std::vector<std::vector<int>> get_distinct_kmer_matrix();
    std::vector<std::vector<int>> get_total_kmer_matrix();
    std::vector<std::vector<int>> get_all_total_matrix();
    std::vector<std::vector<int>> get_all_shared_matrix();
    std::vector<std::vector<int>> get_distinct_total_matrix();
    std::vector<std::vector<int>> get_distinct_shared_matrix();

  private:
    GraphWrapper& graph_ptr;
    std::vector<std::vector<int>> distinct_kmer_matrix;
    std::vector<std::vector<int>> total_kmer_matrix;

    std::vector<std::vector<int>> all_shared_matrix;
    std::vector<std::vector<int>> all_total_matrix;

    std::vector<std::vector<int>> distinct_shared;
    std::vector<std::vector<int>> distinct_total;

    std::vector<std::vector<int>> phenotypes;

    int core_threshold;
    int unique_threshold;
    int genome_number;
    std::map<std::string, std::vector<int>> kmer_map;
};
