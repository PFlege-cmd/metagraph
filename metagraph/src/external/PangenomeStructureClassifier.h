//
// Created by Patrick Flege on 29/09/2026.
//

#ifndef METAGRAPH_PANGENOMESTRUCTURECLASSIFIER_H
#define METAGRAPH_PANGENOMESTRUCTURECLASSIFIER_H
#include "PangenomeStructureAdapter.h"
#include "cli/GraphWrapper.hpp"


class PangenomeStructureClassifier {
    //TODO: Make frequencies union-type, and part of the PangenomeCollector. Flatten std::vector<std::vector<int>> to std::vector<int>
    public:
    typedef PangenomeStructureAdapter::PangenomeStructureCollector StructureCollector;

        PangenomeStructureClassifier(GraphWrapper &graph_wrapper, int number_of_genomes, int loops)
    :
    m_number_of_genomes(number_of_genomes),
    m_loops(loops) {
            std::cout << "PangenomeStructureClassifier constructed" << std::endl;
            std::cout << "Number of genomes: " << m_number_of_genomes << std::endl;
            std::cout << "Loops: " << m_loops << std::endl;
        };

    /**@brief Calculates random combination of collection of genome indices.
     * Example: [0, 1, 2, 3] -> [2, 0, 1, 3]
     * Uses Marsenne-twister engine and uniform distribution to get pseudo-random combinations.
     * @param generator         Marsenne-twister engine for random-number generation.
     * @param distributions     Random uniform distributions, for numbers in the range 'genome_number' down to 2.
     * @return Random combination of indices in the range of the size of the number of genomes.
     */
std::vector<int> generate_genome_combinations(std::mt19937 &generator, std::vector<std::uniform_int_distribution<int>> &distributions);

    /** @brief Generates pseudo-random genome-id combinations for a number of loops (defined at construction).
     *
     * @return vector of 'loop'-pseudorandom combinations of genome-ids.
     */
std::vector<std::vector<int>> generate_genome_combinations_for_loop() const;





    /**
         * @brief For a given loop-iteration and genome (index) combination, calculate the core-, accessory- and unique-kmers of the collector.
         * @param collector Collector for distinct-kmer frequency-categories and size of (sub)-pangenome, contains bit-valued frequencies.
         * @param loop Current loop-iteration
         * @param growing_combinations current subset of powerset of all possible combinations.
    */
void integrate_frequencies_for_combination(
        StructureCollector& collector,
            int loop,
            const std::vector<int>& growing_combinations) const;

/**@brief Loops through a set of random combinations of genome-ids of the pangenome, iteratively classifying it.
 *
 * @param collector Collector for core-, accessory-, and unique kmers
 * @param combinations Combinations of genome-ids.
 * @return Updated collector.
 */
PangenomeStructureAdapter::PangenomeStructureCollector& loop_for_increasing_pangenome_size(StructureCollector& collector, std::vector<std::vector<int>>& combinations) const;

    private:
        int m_number_of_genomes{0};
        int m_loops{0};
};




#endif // METAGRAPH_PANGENOMESTRUCTURECLASSIFIER_H
