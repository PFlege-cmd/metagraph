//
// Created by Patrick Flege on 29/09/2026.
//

#include <random>

#include "PangenomeStructureClassifier.h"

std::vector<int> PangenomeStructureClassifier::generate_genome_combinations(std::mt19937 &generator, std::vector<std::uniform_int_distribution<int>> &distributions) {

    std::vector<int> genome_ids(m_number_of_genomes, 0);
    std::iota(genome_ids.begin(), genome_ids.end(), 0);



    std::vector<int> genome_id_combinations{};
    for (int i = m_number_of_genomes - 1; i >= 1; i--) {
        int next_idx{distributions[i](generator)};

        int next_genome_number{genome_ids[next_idx]};
        genome_id_combinations.push_back(next_genome_number);

        std::vector<int> updated_idxs{};
        std::remove_copy_if(genome_ids.begin(), genome_ids.end(), std::back_inserter(updated_idxs), [next_genome_number](int gen_idx) {
            return gen_idx == next_genome_number;
        });
        genome_ids = updated_idxs;
    }

    return genome_id_combinations;
}


std::vector<std::vector<int>> PangenomeStructureClassifier::generate_genome_combinations_for_loop() {
    std::mt19937 generator{std::random_device{}()};
    std::vector<std::uniform_int_distribution<int>> distributions(m_number_of_genomes);
    for (size_t i = 0; i < static_cast<size_t>(m_number_of_genomes); i++) {
        distributions[i] = std::uniform_int_distribution<int>{0, static_cast<int>(i)};
    }

    std::vector<std::vector<int>> genome_id_combinations(m_loops);

    for (size_t i = 0; i < static_cast<size_t>(m_loops); i++) {
        genome_id_combinations[i] = generate_genome_combinations(generator, distributions);
    }
    return genome_id_combinations;
};
