//
// Created by Patrick Flege on 16/09/2026.
//

#ifndef PANGENOMESTRUCTUREADAPTER_H
#define PANGENOMESTRUCTUREADAPTER_H
#include <vector>

/**
 * @brief Adaptor for connecting to C++-code from PanTools for 'pangenome-structure' command
 *
 * Carries the different list-loop constructs which will be evaluated in PanTools' Java-codebase.
 */
class PangenomeStructureAdapter {
public:
    struct PangenomeStructureCollector {
        std::vector < std::vector<int> > core_loops{std::vector<std::vector<int>>(100)};
        std::vector < std::vector<int> > accessory_loops{std::vector<std::vector<int>>(100)};
        std::vector < std::vector<int> > unique_loops{std::vector<std::vector<int>>(100)};
    };
};



#endif //PANGENOMESTRUCTUREADAPTER_H
