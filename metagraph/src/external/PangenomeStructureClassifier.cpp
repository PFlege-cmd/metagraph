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
    // std::mt19937 generator{std::random_device{}()};
    // std::vector<std::uniform_int_distribution<int>> distributions(m_number_of_genomes);
    // for (size_t i = 0; i < static_cast<size_t>(m_number_of_genomes); i++) {
    //     distributions[i] = std::uniform_int_distribution<int>{0, static_cast<int>(i)};
    // }
    //
    // std::vector<std::vector<int>> genome_id_combinations(m_loops);
    //
    // for (size_t i = 0; i < static_cast<size_t>(m_loops); i++) {
    //     genome_id_combinations[i] = generate_genome_combinations(generator, distributions);
    // }
    // return genome_id_combinations;

    return CombinationUtils::generate_combinations_for_loop(m_loops, m_number_of_genomes);
};

void PangenomeStructureClassifier::integrate_frequencies_for_combination(
    StructureCollector& collector,
        std::vector<std::vector<int>>& frequencies,
        int loop,
        const std::vector<int>& growing_combinations) const {
        for (int freq_idx = 0; freq_idx < static_cast<int>(frequencies.size()); freq_idx++) {

            size_t combinations_idx = growing_combinations.size() - 2;

            auto classified_kmers { KmerUtils::classify_selected_kmers(frequencies[freq_idx],
                                                                       growing_combinations) };
            KmerType classification { std::get<0>(classified_kmers) };
            //KmerUtils::genomes_and_frequencies& genome_freqs = std::get<1>(classified_kmers);

            switch (classification) {
                case KmerType::CORE:
                    // for (std::pair<int, int> genome_freq_ : genome_freqs) {
                    //     //int gen_idx { genome_freq.first };
                    //     //TODO: Fix all these unused frequencies and adapt code (tests etc)
                    //     ///int gen_freq { genome_freq.second };
                    //     //collector.core_loops[loop][gen_idx] = gen_freq;
                    //     BOOST_ASSERT(genome_freq_.second != 0);
                    //     collector.core_loops[loop][combinations_idx] += 1;
                    // }
                    collector.core_loops[loop][combinations_idx] += 1;
                    break;
                case KmerType::ACCESSORY:
                    // for (std::pair<int, int> genome_freq_ : genome_freqs) {
                    //     //int gen_idx { genome_freq.first };
                    //     //int gen_freq { genome_freq.second };
                    //     BOOST_ASSERT(genome_freq_.second != 0);
                    //     collector.accessory_loops[loop][combinations_idx] += 1;
                    // }
                    collector.accessory_loops[loop][combinations_idx] += 1;
                    break;
                case KmerType::UNIQUE:
                    // for (std::pair<int, int> genome_freq_ : genome_freqs) {
                    //     //int gen_idx { genome_freq.first };
                    //     //int gen_freq { genome_freq.second };
                    //     BOOST_ASSERT(genome_freq_.second != 0);
                    //     collector.unique_loops[loop][combinations_idx] += 1;
                    // }
                    collector.unique_loops[loop][combinations_idx] += 1;
                    break;
                case KmerType::NOT_PRESENT:
                    break;
            }
        }
    }

PangenomeStructureClassifier::StructureCollector& PangenomeStructureClassifier::loop_for_increasing_pangenome_size(StructureCollector& collector, std::vector<std::vector<int>> &frequencies, std::vector<std::vector<int>>& combinations) const {
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
                        integrate_frequencies_for_combination(collector, frequencies, loop,
                                                              growing_combinations);
                    }
                }

        }
    }
    return collector;
}
