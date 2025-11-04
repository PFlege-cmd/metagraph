//
// Created by Patrick Flege on 04/11/2025.
//

#ifndef MOCK_ANNOTATED_DBG_H
#define MOCK_ANNOTATED_DBG_H
#include "graph/annotated_dbg.hpp"
#include <gmock/gmock.h>
#include <gtest/gtest_prod.h>
#include <vector>


class mock_annotated_dbg : mtg::graph::AnnotatedDBG {

    public:
    mock_annotated_dbg() : AnnotatedDBG(nullptr, nullptr) {};
    MOCK_METHOD((std::vector<std::array<int, 2>> ),read_mapping_pantools_both_sides, (std::string_view& read, std::string_view& genome, std::vector<int>& sequence_lengths_vector));

};



#endif //MOCK_ANNOTATED_DBG_H
