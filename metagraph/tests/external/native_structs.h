//
// Created by Patrick Flege on 02/09/2026.
//

#ifndef NATIVE_STRUCTS_H
#define NATIVE_STRUCTS_H
#include <iostream>
#include <__ostream/basic_ostream.h>
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
            delete[] core;
            delete[] accessory;
            delete[] unique;
        }
    };

}
#endif //NATIVE_STRUCTS_H
