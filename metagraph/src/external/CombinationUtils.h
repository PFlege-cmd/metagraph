//
// Created by Patrick Flege on 30/09/2026.
//

#ifndef METAGRAPH_COMBINATIONUTILS_H
#define METAGRAPH_COMBINATIONUTILS_H
#include <numeric>
#include <random>
#include <vector>
#include <algorithm>

namespace CombinationUtils {

/**
 *@brief Generates combinations of increasing numbers, from 0 to max_number.
 *
 * @param generator Generator for random-distribution
 * @param distributions Random distribution to generate numbers
 * @param max_number Maximum of possible values
 * @return
 */
inline std::vector<int> generate_combinations(std::mt19937 &generator, std::vector<std::uniform_int_distribution<int>> &distributions, int max_number) {

    std::vector<int> ids(max_number, 0);
    std::iota(ids.begin(), ids.end(), 0); // IDs to draw from in uniform distribution

    std::vector<int> id_combinations{};
    // Updates a set of ID-values (to be taken from successively smaller distributions
    for (int i = max_number - 1; i >= 1; i--) {
        int next_idx{distributions[i](generator)}; // next index from distributions.

        int next_number{ids[next_idx]};
        id_combinations.push_back(next_number);

        std::vector<int> updated_idxs{};
        std::copy_if(ids.begin(), ids.end(), std::back_inserter(updated_idxs), [next_number](int idx) {
            return idx != next_number;
        });
        // std::remove_copy_if(ids.begin(), ids.end(), std::back_inserter(updated_idxs), [next_number](int idx) {
        //     return idx == next_number;
        // });
        ids = updated_idxs;
    }

    return id_combinations;
}

/**
 * @brief Generates combinations of increasing sequences of numbers for a specified maximum, and number of loops.
 *
 * @param loops Number of different combinations
 * @param max_number Maximum value of sequence-element
 * @return vector of shuffled sequences.
 */
inline std::vector<std::vector<int>> generate_combinations_for_loop(int loops, int max_number) {
    std::mt19937 generator{std::random_device{}()};
    std::vector<std::uniform_int_distribution<int>> distributions(max_number);
    for (size_t i = 0; i < static_cast<size_t>(max_number); i++) {
        distributions[i] = std::uniform_int_distribution<int>{0, static_cast<int>(i)};
    }

    std::vector<std::vector<int>> id_combinations(loops);

    for (size_t i = 0; i < static_cast<size_t>(loops); i++) {
        id_combinations[i] = generate_combinations(generator, distributions, max_number);
    }
    return id_combinations;
};
}


#endif // METAGRAPH_COMBINATIONUTILS_H
