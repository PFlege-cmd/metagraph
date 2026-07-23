//
// Created by Patrick Flege on 17/07/2026.
//
#include "external/MatrixEntryFactory.h"
#include "mock_matrix_entry_factory.h"
std::unique_ptr<MockMatrixEntryFactory> MockMatrixEntryFactory::instance = nullptr;

MockMatrixEntryFactory::MockMatrixEntryFactory(const std::string& genome,
                                               int no_seqs,
                                               int* seqs,
                                               GraphWrapper& wrapper)
    : wrapper(wrapper) {
    this->seqs = std::vector<int>(no_seqs);
    this->genome = genome;
    for (int i = 0; i < no_seqs; i++) {
        this->seqs[i] = seqs[i];
    }
    this->wrapper = wrapper;
};

std::unique_ptr<MockMatrixEntryFactory>& MockMatrixEntryFactory::getInstance(std::string genome, int no_seqs, int* seqs, GraphWrapper& wrapper) {
    if (!MockMatrixEntryFactory::instance) {
        MockMatrixEntryFactory::instance = std::make_unique<MockMatrixEntryFactory>(genome, no_seqs, seqs, wrapper);
    }
    return MockMatrixEntryFactory::instance;
}




