//
// Created by Patrick Flege on 15/07/2026.
//

#ifndef MOCK_MATRIX_ENTRY_H
#define MOCK_MATRIX_ENTRY_H
#include "external/MatrixEntry.h"

#include "gmock/gmock-function-mocker.h"
class MockMatrixEntry : public MatrixEntry {
    public:
        explicit MockMatrixEntry(std::string_view& genome, std::string_view& kmer, std::vector<int>& sequence_lengths, GraphWrapper& graph);
        explicit MockMatrixEntry(MockMatrixEntry& other);


        MOCK_METHOD(std::string_view&, get_reverse_kmer, (), (override));
        MOCK_METHOD((std::vector<std::array<int, 2>>), locate, (), (override));
        MOCK_METHOD((std::vector<std::array<int, 2>>), locate_reverse, (), (override));
        MOCK_METHOD((std::vector<std::array<int, 2>>), locate_both_kmers, (), (override));
};
#endif //MOCK_MATRIX_ENTRY_H
