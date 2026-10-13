//
// Created by Patrick Flege on 13/10/2026.
//

#include "external/KmerIDCollector.h"


#include "gtest/gtest.h"
TEST(testKmerIDCollector, testMergeToSource) {

    std::vector<unsigned long long> core_ids{1, 2, 3};
    std::vector<std::vector<unsigned long long>> accessory_ids{{4, 5},{7,8}};
    std::vector<std::vector<unsigned long long>> unique_ids{{9, 10},{11, 12}, {17}, {}};
    KmerIDCollector origin{std::move(core_ids),std::move(accessory_ids), std::move(unique_ids)};

    std::vector<unsigned long long> local_core_ids{111, 22};
    std::vector<std::vector<unsigned long long>> local_accessory_ids{{19, 21},{}};
    std::vector<std::vector<unsigned long long>> local_unique_ids{{9, 10},{11, 12}, {},{13}};
    KmerIDCollector delta{std::move(local_core_ids), std::move(local_accessory_ids), std::move(local_unique_ids)};
    origin.merge(delta);

    ASSERT_EQ(origin.core_ids.size(), 5);

    ASSERT_EQ(origin.accessory_ids[0].size(), 4);
    ASSERT_EQ(origin.accessory_ids[1].size(), 2);

    ASSERT_EQ(origin.unique_ids[0].size(), 4);
    ASSERT_EQ(origin.unique_ids[1].size(), 4);
    ASSERT_EQ(origin.unique_ids[2].size(), 1);
    ASSERT_EQ(origin.unique_ids[3].size(), 1);


    //TODO: Figure out what the moved-on state should be.
    //
    // ASSERT_EQ(delta.accessory_ids[0].size(), 0);
    // ASSERT_EQ(delta.accessory_ids[1].size(), 0);
    //
    // ASSERT_EQ(delta.unique_ids[0].size(), 0);
    // ASSERT_EQ(delta.unique_ids[1].size(), 0);
    // ASSERT_EQ(delta.unique_ids[2].size(), 0);
    // ASSERT_EQ(delta.unique_ids[3].size(), 0);

}

TEST(testKmerIDCollector, testInitialize) {
    auto coll{KmerIDCollector::initialize_kmer_collector<4>()};
    ASSERT_EQ(coll->accessory_ids.size(), 2);
    ASSERT_EQ(coll->unique_ids.size(), 4);
}