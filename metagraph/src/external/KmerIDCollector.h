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
    KmerIDCollector* merge( KmerIDCollector& other) {
        core_ids.insert(core_ids.end(), std::make_move_iterator(other.core_ids.begin()), std::make_move_iterator(other.core_ids.end()));

        size_t accessory_size{other.accessory_ids.size()};
        size_t unique_size{other.unique_ids.size()};

        for (size_t i = 0; i < accessory_size; ++i) {
            //accessory_ids[i].insert(accessory_ids[i].end(), std::make_move_iterator(other.accessory_ids[i].begin()), std::make_move_iterator(other.accessory_ids[i].end()));
            std::move(other.accessory_ids[i].begin(), other.accessory_ids[i].end(), std::back_inserter(accessory_ids[i]));
        }

        for (size_t i = 0; i < unique_size; ++i) {
            unique_ids[i].insert(unique_ids[i].end(), std::make_move_iterator(other.unique_ids[i].begin()), std::make_move_iterator(other.unique_ids[i].end()));
        }
        return this;
    }

    template <int GENOMES>
    static KmerIDCollector* initialize_kmer_collector() {
        std::vector<unsigned long long> core;
        std::vector<std::vector<unsigned long long>> acc (GENOMES - 2, std::vector<unsigned long long>{});
        std::vector<std::vector<unsigned long long>> uniq (GENOMES, std::vector<unsigned long long>{});

        return new KmerIDCollector{ .core_ids = std::move(core), .accessory_ids = std::move(acc), .unique_ids = std::move(uniq) };
    }
};
#endif // METAGRAPH_KMERIDCOLLECTOR_H
