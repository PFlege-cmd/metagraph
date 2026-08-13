//
// Created by Patrick Flege on 01/09/2026.
//

#ifndef READMAPPER_H
#define READMAPPER_H

#include "../../tests/external/native_structs.h"


#include <map>
#include <vector>

#include "cli/GraphWrapper.hpp"
#include "common/seq_tools/reverse_complement.hpp"

class ReadMapper {
public:
    using sequence_and_position = std::array<int, 2>;
    ReadMapper(ReadMapper &) = delete;
    ReadMapper(std::vector<std::vector<unsigned long long>>& sequence_offsets,
        std::vector<std::vector<unsigned long long>>& sequence_lengths,
        std::map <std::string_view, std::vector<sequence_and_position>>& genome_hits_map,
        std::map<std::string_view, size_t>& genome_name_to_index_map,
        GraphWrapper& wrapper, const std::string_view& read)
    : m_sequence_offsets{std::move(sequence_offsets)}, m_sequence_lengths{std::move(sequence_lengths)}, m_genome_hits_map {std::move(genome_hits_map)}, m_genome_name_to_index_map{genome_name_to_index_map}, m_wrapper {wrapper}, m_read {read}, m_reverse_read{m_read} {
        //std::string reverse_read = std::string(m_read);
        reverse_complement(m_reverse_read); // Prevents backing string from disappearing
        std::cout << "Reverse read is: " << m_reverse_read << std::endl;
        //m_reverse_read = reverse_read;
    }
    explicit ReadMapper(GraphWrapper& wrapper) : m_wrapper(wrapper) {};

    std::map<std::string_view, std::vector<sequence_and_position>> &get_genome_hits_map() {
        return m_genome_hits_map;
    }

    void set_read(std::string_view read) {
        m_read = read;
        m_reverse_read = read;
        reverse_complement(m_reverse_read);
    }

    void initialize_all_fields(int* seqLengths, int* seqsPerGenome);
    void read_mapping_both_sides();
    void read_mapping_forward();
    void read_mapping_backward();
    void calculate_sequence_location(unsigned long coord, const std::vector<unsigned long long>& sequence_offsets, const std::vector<unsigned long long>& sequence_lengths, std::array<int, 2>& position_and_location);
    void calculate_sequence_location_reverse(unsigned long coord,const std::vector<unsigned long long>& sequence_offsets, const std::vector<unsigned long long>& sequence_lengths, std::array<int, 2>& position_and_location);
    NativeStructures::HitsPerGenomes* transform_to_native_output();

private:
    std::vector<std::vector<unsigned long long>> m_sequence_offsets{};
    std::vector<std::vector<unsigned long long>> m_sequence_lengths{};
    std::map <std::string_view, std::vector<sequence_and_position>> m_genome_hits_map{};
    std::map<std::string_view, size_t> m_genome_name_to_index_map{};
    GraphWrapper& m_wrapper;
    std::string_view m_read{};
    std::string m_reverse_read{};
    int m_kmer_dist{15};
};
#endif //READMAPPER_H
