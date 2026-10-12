//
// Created by Patrick Flege on 29/09/2026.
//

#include <random>

#include "PangenomeStructureClassifier.h"

#include "KmerClassifier.hpp"
#include "KmerUtils.h"
#include "CombinationUtils.h"

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


std::vector<std::vector<int>> PangenomeStructureClassifier::generate_genome_combinations_for_loop() const {
    return CombinationUtils::generate_combinations_for_loop(m_loops, m_number_of_genomes);
};

void PangenomeStructureClassifier::integrate_frequencies_for_combination(
    StructureCollector& collector,
        int loop,
        const std::vector<int>& growing_combinations) const {
        for (int freq_idx = 0; freq_idx < static_cast<int>(collector.number_nodes); freq_idx++) {

            size_t combinations_idx = growing_combinations.size() - 2;

            auto classification = KmerUtils::presence_absence_kmers(freq_idx, this->m_number_of_genomes, collector.freqs, growing_combinations) ;
            switch (classification) {
                case KmerType::CORE:
                    collector.core_loops[loop][combinations_idx] += 1;
                    break;
                case KmerType::ACCESSORY:
                    collector.accessory_loops[loop][combinations_idx] += 1;
                    break;
                case KmerType::UNIQUE:
                    collector.unique_loops[loop][combinations_idx] += 1;
                    break;
                case KmerType::NOT_PRESENT:
                    break;
            }
        }
    }

PangenomeStructureClassifier::StructureCollector& PangenomeStructureClassifier::loop_for_increasing_pangenome_size(StructureCollector& collector, std::vector<std::vector<int>>& combinations) const {
    // Outermost loop: Loop through each loop-iteration
    for (int loop = 0; loop < m_loops; loop++) {

        int first_entry_id{0};
        // Get the first combination with 2 entries (only one is meaningless)
        std::vector<int>& cur_combinations = combinations[loop];
        std::vector<int> growing_combinations{};
        growing_combinations.push_back(cur_combinations[first_entry_id]);

        #pragma omp parallel
        {
            #pragma omp for schedule(static, 3)
            for (int idx = 1; idx < static_cast<int>(cur_combinations.size()); idx++) {
                    // Innermost loop: Go through all frequencies, update the core-, accessory-, and unique entries for each loop-combination join
                #pragma omp critical
                {
                        growing_combinations.push_back(cur_combinations[idx]);
                        integrate_frequencies_for_combination(collector, loop,
                                                              growing_combinations);
                }
            }
        }
    }
    return collector;
}
