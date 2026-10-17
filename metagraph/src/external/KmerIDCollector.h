//
// Created by Patrick Flege on 13/10/2026.
//

#ifndef METAGRAPH_KMERIDCOLLECTOR_H
#define METAGRAPH_KMERIDCOLLECTOR_H
#include <vector>


/**
 * @brief Collector for Kmer-ids when classifying the kmer-profile of a pangenome.
 *
 * Stores core-ids, accessory-ids (ordered by increasing number of shared genomes), and unique (ordered by the genome-id where they are found).
 *
 * Example: Pangenome has four genomes:
* Core is simply stored in 'core_ids'.
* Accessory is assigned as follows:
*   - Occurring in 2 genomes: First vector of accessory_ids,
*   - Occurring in 3 genomes: Second vector of accessory_ids.
 * Unique are partitioned via the id of their origin genome.
 *  - Only in genome 0 (first genome) : assigned to unique_ids[0]
 *  - Only in genome 1 (second genome) : assigned to unique_ids[1]
 */
struct KmerIDCollector {
    std::vector<unsigned long long> core_ids;
    std::vector<std::vector<unsigned long long> > accessory_ids;
    std::vector<std::vector<unsigned long long> > unique_ids;

    /**
     * @brief Merges the entries of the other, local kmer-collector into the current one. Uses move-semantics.
     * @param other (Local) KmerIDCollector from which ids are drawn.
     * @return Updated (this) collector.
     */
    KmerIDCollector& merge(const KmerIDCollector& other) {
        core_ids.insert(core_ids.end(), std::make_move_iterator(other.core_ids.begin()), std::make_move_iterator(other.core_ids.end()));

        std::size_t accessory_size{other.accessory_ids.size()};
        std::size_t unique_size{other.unique_ids.size()};

        for (std::size_t i = 0; i < accessory_size; ++i) {
            std::move(other.accessory_ids[i].begin(), other.accessory_ids[i].end(), std::back_inserter(accessory_ids[i]));
        }

        for (std::size_t i = 0; i < unique_size; ++i) {
            std::move(other.unique_ids[i].begin(), other.unique_ids[i].end(), std::back_inserter(unique_ids[i]));
        }
        return *this;
    }

    static KmerIDCollector* initialize_kmer_collector(int genomes) {
        std::vector<unsigned long long> core;
        std::vector<std::vector<unsigned long long>> acc (genomes - 2);
        std::vector<std::vector<unsigned long long>> uniq (genomes);

        auto temp = new KmerIDCollector{ .core_ids = std::vector<unsigned long long>{},
            .accessory_ids = std::vector<std::vector<unsigned long long>> (genomes-2),
            .unique_ids = std::vector<std::vector<unsigned long long>> (genomes) };

        return temp;
    }
};
#endif // METAGRAPH_KMERIDCOLLECTOR_H
