//
// Created by Patrick Flege on 17/07/2026.
//

#ifndef GENOMECOORDINATEFACTORY_H
#define GENOMECOORDINATEFACTORY_H
#include "GenomeCoordinate.h"
#include "MatrixEntryFactory.h"


class GenomeCoordinateFactory : public MatrixEntryFactory {
protected:
    static std::unique_ptr<MatrixEntryFactory> instance;
public:
    GenomeCoordinateFactory(std::string genome, int no_seqs, int* seqs, GraphWrapper& wrapper);
  std::vector<std::unique_ptr<MatrixEntry>> createKmerEntries(const char** kmers,
                                                                 int numberOfKmers) override;
    static std::unique_ptr<MatrixEntryFactory>& getInstance(std::string genome, int no_seqs, int* seqs, GraphWrapper& wrapper);


private:
    std::string_view genome;
    std::vector<int> seqs;
    GraphWrapper& wrapper;
};



#endif //GENOMECOORDINATEFACTORY_H
