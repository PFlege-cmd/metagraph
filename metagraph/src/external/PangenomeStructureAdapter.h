//
// Created by Patrick Flege on 16/09/2026.
//

#ifndef PANGENOMESTRUCTUREADAPTER_H
#define PANGENOMESTRUCTUREADAPTER_H
#include "native_structs.h"
#include <vector>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/iostreams/filter/zstd.hpp>

/**
 * @brief Adaptor for connecting to C++-code from PanTools for 'pangenome-structure' command
 *
 * Carries the different list-loop constructs which will be evaluated in PanTools' Java-codebase, as well as frequencies.
 */
class PangenomeStructureAdapter {
public:
    struct PangenomeStructureCollector {
        std::vector < std::vector<int> > core_loops{std::vector<std::vector<int>>(100)};
        std::vector < std::vector<int> > accessory_loops{std::vector<std::vector<int>>(100)};
        std::vector < std::vector<int> > unique_loops{std::vector<std::vector<int>>(100)};
        boost::dynamic_bitset<> freqs
                = boost::dynamic_bitset<>(boost::dynamic_bitset<>(100));
        size_t number_nodes{1};
    };

    /**
     *
     * @param collector Collector for core-, accessory-, and unique entries of kmer-classifi cation combinations.
     * @return Pangenomestructure-matrix with dynamic arrays instead of vectors for export via C-class
     */
    static NativeStructures::PangenomeStructureKmerMatrix*
    export_matrix(const PangenomeStructureCollector& collector) {
        //int genomes{ static_cast<int>(collector.core_loops[0].size()) };
        int combinations{static_cast<int>(collector.core_loops[0].size())};
        int loops{ static_cast<int>(collector.core_loops.size()) };

        long * core_entries = new (std::nothrow) long[combinations*loops]();
        long * accessory_entries = new (std::nothrow) long[combinations*loops]();
        long * unique_entries = new (std::nothrow) long[combinations*loops]();
        for (int loop = 0; loop < loops; loop++) {
            for (int cmb_idx = 0; cmb_idx < combinations; cmb_idx++) {
                core_entries[loop*combinations + cmb_idx] = collector.core_loops[loop][cmb_idx];
                accessory_entries[loop*combinations + cmb_idx] = collector.accessory_loops[loop][cmb_idx];
                unique_entries[loop*combinations + cmb_idx] = collector.unique_loops[loop][cmb_idx];
            }
        }
        return new NativeStructures::PangenomeStructureKmerMatrix{core_entries, accessory_entries, unique_entries};
    }
};



#endif //PANGENOMESTRUCTUREADAPTER_H
