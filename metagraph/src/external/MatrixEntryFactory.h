//
// Created by Patrick Flege on 17/07/2026.
//

#ifndef MATRIXENTRYFACTORY_H
#define MATRIXENTRYFACTORY_H
#include <algorithm>
#include "MatrixEntry.h"

class MatrixEntryFactory {

public:
        MatrixEntryFactory();
        virtual ~MatrixEntryFactory() = default;
        virtual std::vector<std::unique_ptr<MatrixEntry>>
        createKmerEntries(const char** kmers, int numberOfKmers)
                = 0 ;
};



#endif //KMERENTRYFACTORY_H
