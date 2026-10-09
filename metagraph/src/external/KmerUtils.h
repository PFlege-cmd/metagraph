//
// Created by Patrick Flege on 29/09/2026.
//

#ifndef METAGRAPH_KMERUTILS_H
#define METAGRAPH_KMERUTILS_H
#include "KmerType.h"


#include <tuple>
#include <utility>
#include <vector>
namespace KmerUtils {

/**
 * Alias for a vector of a kmer-type and frequency for this.
 */
typedef std::vector<std::pair<int, int>> genomes_and_frequencies;

    /**
     * Method to classify kmers based on their frequency in a selected set of genomes.
     * @param kmer_counts       Counts of kmers across all genomes
     * @param genome_indices    Indices of genomes for which we wish to classify genomes.
     * @return
     */
    static std::tuple<KmerType, genomes_and_frequencies> classify_selected_kmers(const std::vector<int>& kmer_counts, const std::vector<int>& genome_indices) {
        KmerType type;
        auto number_of_occuring_genomes{0};
        genomes_and_frequencies genome_frequencies({});
        int unique_threshold{1};
        int core_threshold{static_cast<int>(genome_indices.size())};

        for (int g_idx: genome_indices) {
            if (kmer_counts[g_idx] > 0 ) {
                number_of_occuring_genomes++;
                auto occurrence_pair = std::make_pair(g_idx, kmer_counts[g_idx]);
                genome_frequencies.push_back(occurrence_pair);
            }
        }

        if (number_of_occuring_genomes == unique_threshold) {
            //std::cout << "Working" << std::endl;
            type = KmerType::UNIQUE;
        }
        else if (number_of_occuring_genomes == core_threshold) {
            //std::cout << "Not unique" << std::endl;
            type = KmerType::CORE;
        }
        else  if (number_of_occuring_genomes > 0) {
            type = KmerType::ACCESSORY;
        } else {
            type = KmerType::NOT_PRESENT;
        }

        return std::tuple<KmerType, genomes_and_frequencies>(type, genome_frequencies);
    };
}
#endif // METAGRAPH_KMERUTILS_H
