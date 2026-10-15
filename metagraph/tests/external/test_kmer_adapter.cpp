//
// Created by Patrick Flege on 14/10/2026.
//

#include <gtest/gtest.h>

#include "../src/external/KmerClassifierAdapter.h"
TEST(test_kmer_adapter, fill_Kmer_Ids) {
    KmerIDCollector* collector = KmerIDCollector::initialize_kmer_collector(4);
    NativeStructures::KmerIds* kmer_ids = new NativeStructures::KmerIds[7]();

    collector->core_ids = std::vector<unsigned long long>{1, 2, 3};

    collector->accessory_ids[0] = std::vector<unsigned long long>{4, 5, 6};
    collector->accessory_ids[1] = std::vector<unsigned long long>{7, 8};

    collector->unique_ids[0] = std::vector<unsigned long long>{9, 10, 1111, 2222, 3444, 5555, 7777, 888};
    collector->unique_ids[1] = std::vector<unsigned long long>{11, 12, 999, 1010, 23232};
    collector->unique_ids[2] = std::vector<unsigned long long>{13, 14};
    collector->unique_ids[3] = std::vector<unsigned long long>{15};

    kmer_ids = KmerClassifierAdapter::fill_Kmer_Ids(*collector, kmer_ids);

    ASSERT_EQ(kmer_ids[1].number_ids, collector->accessory_ids[0].size());
    ASSERT_EQ(kmer_ids[2].number_ids, collector->accessory_ids[1].size());

    ASSERT_EQ(kmer_ids[3].number_ids, collector->unique_ids[0].size());
    ASSERT_EQ(kmer_ids[4].number_ids, collector->unique_ids[1].size());
    ASSERT_EQ(kmer_ids[5].number_ids, collector->unique_ids[2].size());
    ASSERT_EQ(kmer_ids[6].number_ids, collector->unique_ids[3].size());

    for (size_t core_id = 0; core_id < collector->core_ids.size(); core_id++) {
        ASSERT_EQ(kmer_ids[0].ids[core_id], collector->core_ids[core_id]);
    }

    for (size_t accessory_id = 0; accessory_id < collector->accessory_ids[0].size(); accessory_id++) {
        ASSERT_EQ(kmer_ids[1].ids[accessory_id], collector->accessory_ids[0][accessory_id]);
    }

    for (size_t accessory_id = 0; accessory_id < collector->accessory_ids[1].size(); accessory_id++) {
        ASSERT_EQ(kmer_ids[2].ids[accessory_id], collector->accessory_ids[1][accessory_id]);
    }

    for (size_t unique_id = 0; unique_id < collector->unique_ids[0].size(); unique_id++) {
        ASSERT_EQ(kmer_ids[3].ids[unique_id], collector->unique_ids[0][unique_id]);
    }

    for (size_t unique_id = 0; unique_id < collector->unique_ids[1].size(); unique_id++) {
        ASSERT_EQ(kmer_ids[4].ids[unique_id], collector->unique_ids[1][unique_id]);
    }

    for (size_t unique_id = 0; unique_id < collector->unique_ids[2].size(); unique_id++) {
        ASSERT_EQ(kmer_ids[5].ids[unique_id], collector->unique_ids[2][unique_id]);
    }

    for (size_t unique_id = 0; unique_id < collector->unique_ids[3].size(); unique_id++) {
        ASSERT_EQ(kmer_ids[6].ids[unique_id], collector->unique_ids[3][unique_id]);
    }
}