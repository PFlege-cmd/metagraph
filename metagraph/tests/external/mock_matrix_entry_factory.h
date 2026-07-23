//
// Created by Patrick Flege on 17/07/2026.
//

#ifndef MOCK_MATRIX_ENTRY_FACTORY_H
#define MOCK_MATRIX_ENTRY_FACTORY_H
#include "external/MatrixEntryFactory.h"
#include "gmock/gmock-function-mocker.h"


class MockMatrixEntryFactory : public MatrixEntryFactory {
    protected:
        static std::unique_ptr<MockMatrixEntryFactory> instance;
    public:
        MockMatrixEntryFactory(const std::string& genome,
                               int no_seqs,
                               int* seqs,
                               GraphWrapper& wrapper);
        static std::unique_ptr<MockMatrixEntryFactory>& getInstance(std::string genome, int no_seqs, int* seqs, GraphWrapper& wrapper);

        MOCK_METHOD((std::vector<std::unique_ptr<MatrixEntry>>), createKmerEntries, (const char**, int), (override));
private:
    std::string_view genome;
    std::vector<int> seqs;
    GraphWrapper& wrapper;
};



#endif //MOCK_MATRIX_ENTRY_FACTORY_H
