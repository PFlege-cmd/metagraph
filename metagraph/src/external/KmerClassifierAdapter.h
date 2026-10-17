//
// Created by Patrick Flege on 22/04/2026.
//

#ifndef KMERCLASSIFIERADAPTER_H
#define KMERCLASSIFIERADAPTER_H
#include "KmerClassifier.hpp"
#include "KmerIDCollector.h"


class KmerClassifierAdapter {
    public:
    /**
     * @brief Exports the results of the kmer-classification into a native-matrix, to be used in PanTools.
     * @param kmerClassifier Classifier of kmer-profile, contains all necessary results.
     */
    static NativeStructures::KmerMatrix* export_matrix(KmerClassifier& kmerClassifier) {
        auto total_kmer_matrix = kmerClassifier.get_total_kmer_matrix();
        auto distinct_kmer_matrix = kmerClassifier.get_distinct_kmer_matrix();

        for (size_t i = 0; i < total_kmer_matrix.size(); i++) {
            for (size_t j = 0; j < total_kmer_matrix[i].size(); j++) {
                std::cout << "Kmer entry at :" << i << j << "--" << total_kmer_matrix[i][j] << std::endl;
            }
        }
        int* core_total_array = new int[total_kmer_matrix[0].size()];
        std::move(total_kmer_matrix[0].begin(), total_kmer_matrix[0].end(), core_total_array);

        int* accessory_total_array = new int[total_kmer_matrix[1].size()];
        std::move(total_kmer_matrix[1].begin(), total_kmer_matrix[1].end(), accessory_total_array);

        int* unique_total_array = new int[total_kmer_matrix[2].size()];
        std::move(total_kmer_matrix[2].begin(), total_kmer_matrix[2].end(), unique_total_array);

        int* core_distinct_array = new int[distinct_kmer_matrix[0].size()];
        std::move(distinct_kmer_matrix[0].begin(), distinct_kmer_matrix[0].end(), core_distinct_array);

        int* accessory_distinct_array = new int[distinct_kmer_matrix[1].size()];
        std::move(distinct_kmer_matrix[1].begin(), distinct_kmer_matrix[1].end(), accessory_distinct_array);

        int* unique_distinct_array = new int[distinct_kmer_matrix[2].size()];
        std::move(distinct_kmer_matrix[2].begin(), distinct_kmer_matrix[2].end(), unique_distinct_array);

        std::vector<std::vector<int>> distinct_shared = kmerClassifier.get_distinct_shared_matrix();
        std::vector<std::vector<int>> distinct_total = kmerClassifier.get_distinct_total_matrix();
        std::vector<std::vector<int>> all_shared = kmerClassifier.get_all_shared_matrix();
        std::vector<std::vector<int>> all_total = kmerClassifier.get_all_total_matrix();

        KmerIDCollector kmerIDCollector{kmerClassifier.get_id_collector()};

        size_t kmer_id_arrays_size{1 + kmerIDCollector.accessory_ids.size() +  kmerIDCollector.unique_ids.size()};

        NativeStructures::KmerIds* kmerIds = new NativeStructures::KmerIds[kmer_id_arrays_size]();

        kmerIds = KmerClassifierAdapter::fill_Kmer_Ids(kmerIDCollector, kmerIds);

        NativeStructures::KmerMatrix *kmerMatrix = new NativeStructures::KmerMatrix(core_total_array,
            accessory_total_array,
            unique_total_array,
            core_distinct_array,
            accessory_distinct_array,
            unique_distinct_array,
            kmerClassifier.flatten(distinct_shared),
            kmerClassifier.flatten(distinct_total),
            kmerClassifier.flatten(all_shared),
            kmerClassifier.flatten(
                        all_total), kmerIds, nullptr, nullptr, nullptr,
                        kmerClassifier.get_total_number_kmers(), kmerClassifier.get_distinct_number_kmers()
            );
        return kmerMatrix;
    };

    static NativeStructures::KmerIds* fill_Kmer_Ids(KmerIDCollector& collector, NativeStructures::KmerIds* kmerIds) {
        long long * core_ids = new long long[collector.core_ids.size()];
        for (size_t i = 0; i < collector.core_ids.size(); i++) {
            core_ids[i] = collector.core_ids[i];
        }

        kmerIds[0] = NativeStructures::KmerIds{static_cast<int>(collector.core_ids.size()), core_ids };
        size_t current_position_kmer_ids{0};

        for (size_t accessory_arrays = 0; accessory_arrays < collector.accessory_ids.size(); accessory_arrays++) {
            current_position_kmer_ids++;
            long long* accessory_ids = new long long[collector.accessory_ids[accessory_arrays].size()];
            kmerIds[current_position_kmer_ids] = NativeStructures::KmerIds{static_cast<int>(collector.accessory_ids[accessory_arrays].size()) , accessory_ids};

            for (size_t accessory_size = 0; accessory_size < collector.accessory_ids[accessory_arrays].size(); accessory_size++) {
                accessory_ids[accessory_size] = collector.accessory_ids[accessory_arrays][accessory_size];
            }
        }

        for (size_t unique_arrays = 0; unique_arrays < collector.unique_ids.size(); unique_arrays++) {
            current_position_kmer_ids++;
            long long* unique_ids = new long long[collector.unique_ids[unique_arrays].size()];
            kmerIds[current_position_kmer_ids] = NativeStructures::KmerIds{static_cast<int>(collector.unique_ids[unique_arrays].size()), unique_ids};
            for (size_t unique_size = 0; unique_size < collector.unique_ids[unique_arrays].size(); unique_size++) {
                unique_ids[unique_size] = collector.unique_ids[unique_arrays][unique_size];
            }
        }

        return kmerIds;
    }
};



#endif //KMERCLASSIFIERADAPTER_H
