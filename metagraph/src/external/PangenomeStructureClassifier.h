//
// Created by Patrick Flege on 29/09/2026.
//

#ifndef METAGRAPH_PANGENOMESTRUCTURECLASSIFIER_H
#define METAGRAPH_PANGENOMESTRUCTURECLASSIFIER_H
#include "cli/GraphWrapper.hpp"


class PangenomeStructureClassifier {
    public:
        PangenomeStructureClassifier(GraphWrapper &graph_wrapper, int number_of_genomes, int loops)
    : m_graph_wrapper_(graph_wrapper),
    m_number_of_genomes(number_of_genomes),
    m_loops(loops) {
            std::cout << "PangenomeStructureClassifier constructed" << std::endl;
            std::cout << "Number of genomes: " << m_number_of_genomes << std::endl;
            std::cout << "Loops: " << m_loops << std::endl;
            std::cout << "Size graph: " << sizeof(m_graph_wrapper_) << std::endl;
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
std::vector<std::vector<int>> generate_genome_combinations_for_loop();

    private:
        GraphWrapper &m_graph_wrapper_;
        int m_number_of_genomes{0};
        int m_loops{0};
};


#endif // METAGRAPH_PANGENOMESTRUCTURECLASSIFIER_H
