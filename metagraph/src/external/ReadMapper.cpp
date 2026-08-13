//
// Created by Patrick Flege on 01/09/2026.
//

#include "ReadMapper.h"

#include "asio/detail/string_view.hpp"
void ReadMapper::read_mapping_both_sides() {

};
void ReadMapper::read_mapping_forward() {
    using Label = std::string;

    int read_length = static_cast<int>(m_read.length());
    std::vector<int> candidate_coords;
    std::vector<std::array<int, 2>> kmer_positions;

    int step_size = read_length / m_kmer_dist;
    step_size = (step_size == 0 ? 1 : step_size);

    int k = static_cast<int>(m_wrapper.get_K());
    for (int position = 0; position < (int) m_read.size(); position += step_size) {
        std::string_view current_kmer = m_read.substr(position, k);
        if ((int) current_kmer.size() < k)
            break;
        //std::cout << "STARTING FORWARD READMAPPING!" << std::endl;
        std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>> coordinates = this->m_wrapper.get_kmer_coordinates(current_kmer);
        for (unsigned long coord_idx = 0; coord_idx < coordinates.size();coord_idx++) {
            std::string_view current_genome = std::get<0>(coordinates[coord_idx]);

            size_t curr_genome_idx{m_genome_name_to_index_map.at(current_genome)};
            //std::cout << "Current genome for these coordinates: " << current_genome << std::endl;
            //std::cout << "Current INDEX for these coordinates: " << curr_genome_idx << std::endl;

            auto coords_for_genome = std::get<2>(coordinates[coord_idx]);
            auto coords_for_genome_size = coords_for_genome[0].size();
            for (unsigned long current_coord_idx = 0; current_coord_idx < coords_for_genome_size; current_coord_idx++) {

                std::array<int,2>  target_and_position = {0, 0};
                calculate_sequence_location(coords_for_genome[0][current_coord_idx],
                                            m_sequence_offsets[curr_genome_idx], m_sequence_lengths[curr_genome_idx], target_and_position);
                int target_sequence = target_and_position[1]; // TODO: Make it inline with pantools
                int pantools_location = target_and_position[0];
                int loc = pantools_location - position;

                // unsigned long long curr_seq_length{curr_genome_idx > 0 ?
                //     m_sequence_offsets[curr_genome_idx][target_sequence] - m_sequence_offsets[curr_genome_idx][target_sequence - 1]:
                //     m_sequence_offsets[curr_genome_idx][target_sequence]};

                unsigned long long curr_seq_length{m_sequence_lengths[curr_genome_idx][target_sequence]};

                if (loc >= 0 && static_cast<unsigned long long>(loc) <= curr_seq_length - static_cast<unsigned long long>(read_length)) {
                    candidate_coords.push_back(loc);
                    m_genome_hits_map.at(current_genome).push_back({target_sequence + 1, static_cast<int>(loc)});
                    //std::cout << "At position: " << position << std::endl;
                    //std::cout << "At genome: " << current_genome << std::endl;

                }
            }
        }
    }
};
void ReadMapper::read_mapping_backward() {
//TODO: Requires me to also give in sequence start and ends. This then will convert it totally to PanTools convertible input.
    //TODO: For now, assume that a sequence is a genome.
    using Label = std::string;

    int read_length = m_read.length();
    std::vector<int> candidate_coords;
    std::vector<std::array<int, 2>> kmer_positions;

    int step_size = read_length / m_kmer_dist;
    step_size = (step_size == 0 ? 1 : step_size);
    //std::cout << "Step size is: " << step_size << std::endl;
    //std::cout << "Complete reverse read:" << m_reverse_read << std::endl;
    std::string_view reverse_read_view = m_reverse_read;

    int k = static_cast<int>(m_wrapper.get_K());;
    for (int position = static_cast<int>((m_read.size())) - k; position > 0; position -= step_size) { // Error was here:
        std::string_view current_kmer = reverse_read_view.substr(position, k);
        //std::cout << "Position: " << position << std::endl;
        //std::cout << "Current kmer is: " << current_kmer << std::endl;
        ////std::cout << current_kmer <<  ": is current kmer!" << std::endl;
        if ((int) current_kmer.size() < k) {
            ////std::cout << "Break for KMER: " <<  current_kmer << std::endl;
            break;
        }

        ////std::cout << current_kmer << std::endl;
        std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>> coordinates = this->m_wrapper.get_kmer_coordinates(current_kmer);
        ////std::cout << "Coordinate size: "<< coordinates.size() << std::endl;
        //std::cout << "BACKWARDS " << std::endl;
        for (unsigned long coord_idx = 0; coord_idx < coordinates.size();coord_idx++) {
            std::string_view current_genome = std::get<0>(coordinates[coord_idx]);

            size_t curr_genome_idx{m_genome_name_to_index_map[current_genome]};
            //std::cout << "Current genome for these coordinates: " << current_genome << std::endl;

            auto coords_for_genome = std::get<2>(coordinates[coord_idx]);
            auto coords_for_genome_size = coords_for_genome[0].size();
            for (unsigned long current_coord_idx = 0; current_coord_idx < coords_for_genome_size; current_coord_idx++) {

                std::array<int,2>  target_and_position = {0, 0};
                calculate_sequence_location_reverse(coords_for_genome[0][current_coord_idx],m_sequence_offsets[curr_genome_idx], m_sequence_lengths[curr_genome_idx],
                                             target_and_position);
                int target_sequence = target_and_position[1]; // TODO: Make it inline with pantools
                int pantools_location = target_and_position[0];
                int loc = pantools_location - position +1; //TODO: THis is very basic. I need to read up on readmapping, and testing in C++. AND clean code

                unsigned long long curr_seq_length{ m_sequence_lengths[curr_genome_idx][target_sequence]};


                if (loc >= 0 && static_cast<unsigned long long>(loc) <= curr_seq_length - static_cast<unsigned long long> (read_length)) {
                    candidate_coords.push_back(loc);
                    //std::cout << "At (reversed) position: " << position << std::endl;
                    m_genome_hits_map[current_genome].push_back({target_sequence + 1, -loc });\
                    //std::cout << "At (reversed) genome: " << current_genome << std::endl;

                } else {
                }
            }
        }
    }
};

void ReadMapper::calculate_sequence_location(unsigned long coord, const std::vector<unsigned long long>& sequence_offsets, const std::vector<unsigned long long>& sequence_lengths, std::array<int, 2>& position_and_location) {
    unsigned long current_length = 0;
    int k = m_wrapper.get_K();
    int target_sequence = -1;
    for (int i = 0; i < static_cast<int>(sequence_lengths.size()); i++) {
        //std::cout << "Seq offsets sizes:" << std::endl;
        //std::cout << sequence_offsets.size() << std::endl;
        //int seq_length = sequence_lengths[i] - k + 1;
        int seq_offset = static_cast<int>(sequence_offsets[i]);
        int seq_length = static_cast<int>(sequence_lengths[i]);
        //std::cout << "current offset " << seq_offset << std::endl;
        //std::cout << "current length " << seq_length << std::endl;
        //std::cout << "Current coordinate: " << coord << std::endl;
        current_length = seq_offset;
        if (coord >= current_length && coord < static_cast<unsigned long>(seq_offset + seq_length - k + 1)) {
            target_sequence = i;
            //std::cout << "target sequence " << target_sequence << std::endl;
            break;
        }
        //current_length = seq_offset; // always keeps the start of the current sequence.
        //current_length += seq_length;
    }

    int pantools_location = coord - current_length;
    position_and_location[0] = pantools_location;
    position_and_location[1] = target_sequence;
    //std::cout << "pantools_location " << pantools_location << std::endl;
    //std::cout << "pantools_sequence " << target_sequence << std::endl;
}


void ReadMapper::calculate_sequence_location_reverse(unsigned long coord,const std::vector<unsigned long long>& sequence_offsets, const std::vector<unsigned long long>& sequence_lengths, std::array<int, 2>& position_and_location) {
    unsigned long long current_length = 0;
    int k = static_cast<int>(m_wrapper.get_K());
    int target_sequence = -1;
    for (int i = 0; i < static_cast<int>(sequence_lengths.size()); i++) {
        current_length = sequence_offsets[i];
        //std::cout << "current offset " << current_length << std::endl;
        //std::cout << "current coord in mappper: " << coord << std::endl;
        int seq_length = static_cast<int>(sequence_lengths[i]);
        if (coord >= sequence_offsets[i] && coord < current_length + seq_length - k + 1) {
            target_sequence = i;
            break;
        }
        //current_length = ;
    }

    // for (int i = 0; i < (int)sequence_lengths.size(); i++) {
    //     int seq_length = sequence_lengths[i] - k + 1;
    //     if ((int) coord >= current_length && (int) coord < current_length + seq_length) {
    //         target_sequence = i;
    //         break;
    //     }
    //     current_length += seq_length;
    // }

    int pantools_location = static_cast<int>(coord - current_length);
    position_and_location[0] = pantools_location;
    position_and_location[1] = target_sequence;
}


NativeStructures::HitsPerGenomes* ReadMapper::transform_to_native_output() {

    auto hits_per_sequence_array = new NativeStructures::HitsPerSequence[m_genome_name_to_index_map.size()];
    // for (int i = 0; i < static_cast<int>(m_genome_name_to_index_map.size()); i++) {
    //     //std::cout << "Initialize from readmapper... " << i << std::endl;
    //     (&hits_per_sequence_array)[i] = new NativeStructures::HitsPerSequence{};
    // }
    NativeStructures::HitsPerGenomes* hpg =  new NativeStructures::HitsPerGenomes(static_cast<int>(m_genome_name_to_index_map.size()));

    //std::vector<std::array<int, 2>> foundInGenome = graph->read_mapping_pantools_both_sides(read, genome, sequence_lengths_vector);
    for (const auto&[genome_name, foundInGenome] : m_genome_hits_map) {
        if (foundInGenome.empty())
            continue;

        int total_hits_per_genome = 0;
        size_t genome_idx{m_genome_name_to_index_map.at(genome_name)};
        size_t no_of_sequences{m_sequence_lengths[genome_idx].size()};
        auto *results = &hits_per_sequence_array[genome_idx];
        ////std::cout << "Genome index: " << genome_idx<<  std::endl;
        ////std::cout << "HitsPerSequence::HitsPerSequence address:" << &results <<std::endl; // As i thought, same object...
        ////std::cout << "Readmapping DONE!" << std::endl;
        //results = new HitsPerSequence[no_of_sequences];
        std::vector<int> tmp_hits;
        std::vector<int> seq_offsets;
        //*results = new HitsPerSequence();
        //int* offsets = new int[no_of_sequences];
        ////std::cout << "Current genome:" << genome_name << std::endl; // stuff WHEN NOTHIONG HERE FIX.
        for (int i = 0; i < static_cast<int>(no_of_sequences); i++) { // Goes through all sequences, and adds results...
            ////std::cout << sequence_lengths[i] << std::endl;
            int counter = 0;
            for (int j = 0; j < (int)foundInGenome.size(); j++) {
                ////std::cout << "SEQUENCE: " << foundInGenome[j][0] << std::endl;
                if (foundInGenome[j][0] == i + 1) { // Checks if the offsets
                    tmp_hits.push_back(foundInGenome[j][1]);
                    ////std::cout << "Pushing back..." << foundInGenome[j][1] << std::endl;
                    counter++;
                    ////std::cout << "Sequence: " << foundInGenome[j][0]
                    //          << ", Position: " << foundInGenome[j][1] << std::endl;
                }
            }
            seq_offsets.push_back(counter);
            total_hits_per_genome += counter;
        }
        hpg->hitFreqPerGenome[genome_idx] = total_hits_per_genome;

        ////std::cout << "first temp hits: " << tmp_hits[0] << std::endl;
        results->hitsPerSequence = new int[tmp_hits.size()];
        results->sequenceOffsets = new int[seq_offsets.size()];
        std::copy(tmp_hits.begin(), tmp_hits.end(), results->hitsPerSequence);
        std::copy(seq_offsets.begin(), seq_offsets.end(), results->sequenceOffsets);
        //results->hitsPerSequence = arraymarray;
       // results->sequenceOffsets = offsets;
    }
    hpg->hitsPerSequence = hits_per_sequence_array;
    return hpg;
};

void ReadMapper::initialize_all_fields(int* seqLengths, int* seqsPerGenome) {
    auto& output = m_wrapper.get_genome_labels();

    std::map<std::string_view, std::vector<sequence_and_position>> genome_map{};
    std::map<std::string_view, size_t> genome_name_to_index_map{};

    size_t curr_index{0};
    for_each(output.begin(), output.end(), [&](const auto& label) {
        //std::cout << "Label is:" << label << std::endl;
        genome_map.emplace(label, std::vector<sequence_and_position>());
        genome_name_to_index_map.emplace(label, curr_index);
        curr_index++;
    });

    size_t size = output.size();
    size_t prev_seqs = 0;
    std::vector<std::vector<unsigned long long>> sequenceOffsets{};
    std::vector<std::vector<unsigned long long>> sequenceLengths{};

    //std::cout << "Number of labels in the pangenome: " <<  size << std::endl;
    size_t k{m_wrapper.get_K()};  // necessary for calculating seq correctly

    for (size_t i = 0; i < size; i++) {
        //std::cout << "Current genome: " << i << std::endl;
        unsigned long long seqOffsets{0};
        std::vector<unsigned long long> seqs_offsets{}; // reset for each.
        std::vector<unsigned long long> lengths{}; // also reset for each.

        //std::cout << "Sequences per genome: "<< seqsPerGenome[i] << std::endl;
        for (size_t j = prev_seqs; j < prev_seqs + seqsPerGenome[i]; j++) {
            seqs_offsets.push_back(seqOffsets);
            lengths.push_back(seqLengths[j]); //updated lengths vector per genome
            //std::cout << "seq lengts" << std::endl;
            //std::cout << "Number of sequences: " << j << std::endl;
            //std::cout << seqLengths[j] << std::endl;
            seqOffsets += seqLengths[j];
            seqOffsets = static_cast<unsigned long long>(seqOffsets - k + 1);
        }
        sequenceOffsets.push_back(seqs_offsets); // total seq-offsets
        sequenceLengths.push_back(lengths); // total sequences-lengths
        prev_seqs += seqsPerGenome[i];
    }
    m_sequence_lengths = sequenceLengths;
    m_sequence_offsets = sequenceOffsets;
    m_genome_name_to_index_map = genome_name_to_index_map;
    m_genome_hits_map = genome_map;
}
