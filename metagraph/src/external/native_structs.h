//
// Created by Patrick Flege on 02/09/2026.
//

#ifndef NATIVE_STRUCTS_H
#define NATIVE_STRUCTS_H
#include <iostream>
namespace NativeStructures {
    struct HitsPerSequence {
        int* hitsPerSequence{nullptr};
        int* sequenceOffsets{nullptr};
        HitsPerSequence(){
            //std::cout << "Hit sPer Sequence initialized" << std::endl;
            hitsPerSequence = new int[1]{-1};
            sequenceOffsets = new int[1]{-1};
        };
        HitsPerSequence(const HitsPerSequence& other) : hitsPerSequence(other.hitsPerSequence), sequenceOffsets(other.sequenceOffsets) {
            //std::cout << "HitsPerSequence::HitsPerSequence COPY CONSTRUCTOR" << std::endl;
        };
        explicit HitsPerSequence(int* arr, int* offsets ) { hitsPerSequence = arr; sequenceOffsets = offsets; };
        ~HitsPerSequence() {
            //std::cout << "HitsPerSequence::HitsPerSequence DESTRUCTOR" << std::endl;
            delete[] hitsPerSequence;
            delete[] sequenceOffsets;
        }
    };


    struct HitsPerGenomes {
        int* hitFreqPerGenome{};
        HitsPerSequence* hitsPerSequence{};
        HitsPerGenomes() = default;
        explicit HitsPerGenomes(int number_genomes) : hitFreqPerGenome(new int[number_genomes]) {
            //std::cout << "HitsPerGenome::HitsPerSequence GENOME INIT CONSTRUCTOR" << std::endl;
            for (int i = 0; i < number_genomes; i++) {
                hitFreqPerGenome[i] = 0;
            }
        }
        explicit HitsPerGenomes(HitsPerSequence* arr) { hitsPerSequence = arr; };
        ~HitsPerGenomes() {
            //std::cout << "HitsPerGenomes::HitsPerGenome DESTRUCTOR" << std::endl;
            delete[] hitFreqPerGenome;
            delete[] hitsPerSequence;
        }
    };

    /**
     *@brief Structure to transfer the kmer-classified pangenome to PanTools.
     *
     * This structure returns kmers classified as core (in all genomes evaluated),
     * accessory (not in all, but found in more than one), and unique (found in only one).
     * Documentation can be found here:
     *
     * @link <a
     *href=https://pantools.readthedocs.io/en/stable/analysis/classification.html#pangenome-structure>Pangenome
     *Structure</a>)
     *
     */
    struct PangenomeStructureKmerMatrix {
        long * core{};
        long * accessory{};
        long * unique{};

        ~PangenomeStructureKmerMatrix() {
            std::cout << "~PangenomeStructureKmerMatrix destructor" << std::endl;
            delete[] core;
            delete[] accessory;
            delete[] unique;
        }

        PangenomeStructureKmerMatrix(long* core, long* accessory, long* unique) {
            this->core = core;
            this->accessory = accessory;
            this->unique = unique;
        }

        PangenomeStructureKmerMatrix(const PangenomeStructureKmerMatrix& other) {
            std::cout << "PangenomeStructureKmerMatrix copy constructor" << std::endl;
            core = other.core;
            accessory = other.accessory;
            unique = other.unique;
        }
    };
    struct KmerIds {
        int number_ids;
        long long* ids;
    };

    struct KmerMatrix {
        int* core_total;
        int* accessory_total;
        int* unique_total;

        int* core_distinct;
        int* accessory_distinct;
        int* unique_distinct;

        int* distinct_shared;
        int* distinct_total;
        int* all_shared;
        int* all_total;

        KmerIds* ids;

        int* shared_phenotypes;
        int* specific_phenotypes;
        int* exclusive_phenotypes;
        explicit KmerMatrix(int* core_total_ptr,
                            int* accessory_total_ptr,
                            int* unique_total_ptr,
                            int* core_distinct_ptr,
                            int* accessory_distinct_ptr,
                            int* unique_distinct_ptr,
                            int* distinct_shared_ptr,
                            int* distinct_total_ptr,
                            int* all_shared_ptr,
                            int* all_total_ptr,
                            KmerIds* ids,
                            int* shared_phenotypes_ptr,
                            int* specific_phenotypes_ptr,
                            int* exclusive_phenotypes_ptr)
            :
    core_total(core_total_ptr),
    accessory_total(accessory_total_ptr),
    unique_total(unique_total_ptr),
    core_distinct(core_distinct_ptr),
    accessory_distinct(accessory_distinct_ptr),
    unique_distinct(unique_distinct_ptr),
    distinct_shared(distinct_shared_ptr),
    distinct_total(distinct_total_ptr),
    all_shared(all_shared_ptr),
    all_total(all_total_ptr),
    ids(ids),
    shared_phenotypes(shared_phenotypes_ptr),
    specific_phenotypes(specific_phenotypes_ptr),
    exclusive_phenotypes(exclusive_phenotypes_ptr)
        {
        };
    };
}
#endif //NATIVE_STRUCTS_H
