//
// Created by Patrick Flege on 26/11/2025.
//

#include "DeBruijnGraphWrapper.h"

 DeBruijnGraphWrapper::DeBruijnGraphWrapper(mtg::graph::AnnotatedDBG &graph) : GraphWrapper(&graph) {}

DeBruijnGraphWrapper::~DeBruijnGraphWrapper() {
 }

long DeBruijnGraphWrapper::retrieveAnchorCoordinates(std::string anchor_sequence,
                                                                long anchor_position,
                                                                std::string genome_name) {
    return 0;
}
uint64_t DeBruijnGraphWrapper::retrieveAnchorId(std::string anchor_sequence) {
    return 0;
}
const char * DeBruijnGraphWrapper::get_sequence_for_coords(std::string genome, unsigned long long start, unsigned long long end){
    return "o";
}

uint64_t DeBruijnGraphWrapper::get_first_node_of_coord_range(uint64_t anchor_index, long long start) {

     uint64_t current_kmer = anchor_index;
     std::vector<uint64_t> outgoing_nodes;

     long long next_coordinate = start + 1;

     while (next_coordinate != start ) {
         bool found_next_coordinate = false;
         bool one_outgoing = this->get_graph()->get_graph().has_single_outgoing(current_kmer);
         bool multi_outgoing = this->get_graph()->get_graph().has_multiple_outgoing(current_kmer);
         std::string one_outgoing_string = (one_outgoing)?"True":"False";
         std::string multiple_outgoing_string = (multi_outgoing)?"True":"False";
         //std::cout << "Has one outgoing? " << one_outgoing_string << std::endl;
         //std::cout << "Has multiple outgoing? " << multiple_outgoing_string << std::endl;
         if (!one_outgoing && !multi_outgoing) {
             std::cout << "Edge stuck!" << std::endl;
             std::cout << next_coordinate << std::endl;
             return -1;
         }

         if (multi_outgoing) {
             std::cout << "Multiple outgoing in Search for kmer: " << this->get_graph()->get_graph().get_node_sequence(current_kmer) << std::endl;
             this->get_graph()->get_graph().adjacent_outgoing_nodes(current_kmer, [&](auto i) {outgoing_nodes.push_back(i);});
             for (unsigned long i = 0; i < outgoing_nodes.size(); ++i) {
                 //std::cout << outgoing_nodes[i] << std::endl;
                 //std::cout << this-> get_graph().get_node_sequence(outgoing_nodes[i]) << std::endl;
                 // create std::vector for this.
                 std::vector<u> edge_node = {outgoing_nodes[i]};
                 std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>> edge_node_coords = this->get_graph()->get_kmer_coordinates(edge_node, num_top_labels, discovery_fraction, presence_fraction);
                 for (unsigned long k = 0; k < edge_node_coords.size(); ++k) {
                     std::string curr_genome = std::get<0>(edge_node_coords[k]);
                     if (std::get<0>(edge_node_coords[k]) != genome ) {
                         continue;
                     }
                     auto x = std::get<2>(edge_node_coords[k])[0];
                     bool found  = std::binary_search(x.begin(), x.end(), next_coordinate);
                     if (!found)
                         continue;

                     //std::cout << "Size of range of coordinates: " << std::get<2>(edge_node_coords[k])[0].size()<< std::endl;
                     auto extracted_coords = std::get<2>(edge_node_coords[k]);
                     auto coordinate_size = extracted_coords[0].size();
                     for (unsigned long j =0; j < coordinate_size; ++j) {
                         auto all_coords = extracted_coords[0];
                         auto outgoing_edge_coords = extracted_coords[0][j];
                         //std::cout << outgoing_edge_coords << std::endl;
                         if (outgoing_edge_coords == next_coordinate) {
                             std::cout << "Found Next kmer: " << std::endl;
                             std::cout << this->get_graph()->get_graph().get_node_sequence(outgoing_nodes[i]) << std::endl;
                             std::cout << "Position: " << extracted_coords[0][j] << std::endl;
                             next_coordinate++;
                             //auto next_char = this->get_graph()->get_graph().get_node_sequence(outgoing_nodes[i]).at(size_kmer - 1);
                             found_next_coordinate = true;
                             current_kmer = outgoing_nodes[i];
                             outgoing_nodes.clear();
                             //starting_kmer.append( 1, next_char);
                             break;
                         }
                     }
                     if (found_next_coordinate)
                         break;
                 }
                 if (found_next_coordinate)
                     break;
             }
         } else if (one_outgoing) {
             std::cout << "SINGLE outgoing in Search for kmer: " << this->get_graph()->get_graph().get_node_sequence(current_kmer) << std::endl;
             this->get_graph()->get_graph().adjacent_outgoing_nodes(current_kmer, [&](auto i) {outgoing_nodes.push_back(i);});
             for (unsigned long i = 0; i < outgoing_nodes.size(); ++i) {
                 //std::cout << outgoing_nodes[i] << std::endl;
                 //std::cout << this-> get_graph().get_node_sequence(outgoing_nodes[i]) << std::endl;
                 // create std::vector for this.
                 std::vector<node_index> edge_node = {outgoing_nodes[i]};
                 std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>> edge_node_coords = this->get_graph()->get_kmer_coordinates(edge_node, num_top_labels, discovery_fraction, presence_fraction);
                 for (unsigned long k = 0; k < edge_node_coords.size(); ++k) {
                     std::string curr_genome = std::get<0>(edge_node_coords[k]);
                     /*if (std::get<0>(edge_node_coords[k]) != genome ) {
                         std::cout << "Not correct genome: " << std::endl;
                         std::cout << curr_genome<< std::endl;
                         continue;
                     }*/ //TODO: Fix the degenerate issue here!
                     auto x = std::get<2>(edge_node_coords[k])[0];
                     //bool found  = binary_search(x.begin(), x.end(), next_coordinate);
                     /*
                     if (!found) {
                         std::cout << "Did not find next coordinate" << std::endl;
                         continue;
                     }*/ // TODO: All not necessary with unitig


                     for (unsigned long j =0; j < std::get<2>(edge_node_coords[k])[0].size(); ++j) {
                         auto next_coor = std::get<2>(edge_node_coords[k])[j][0];
                         std::cout << next_coor << std::endl;

                         //std::binary_search(std::get<2>(edge_node_coords[k])[0])
                         //if (std::get<2>(edge_node_coords[k])[0][j] == next_coordinate) {
                         std::cout << "Found Next kmer: " << std::endl;
                         std::cout << this->get_graph()->get_graph().get_node_sequence(outgoing_nodes[i]) << std::endl;
                         std::cout << "Position: " << std::get<2>(edge_node_coords[k])[0][j] << std::endl;
                         next_coordinate++;
                         //char next_char = this->get_graph()->get_graph().get_node_sequence(outgoing_nodes[i]).at(size_kmer - 1);
                         found_next_coordinate = true;
                         current_kmer = outgoing_nodes[i];
                         outgoing_nodes.clear();
                         //starting_kmer.append( 1, next_char);
                         break;
                         //}
                     }
                     if (found_next_coordinate)
                         break;
                 }
                 if (found_next_coordinate)
                     break;
             }
         }
         return 0;
     }
 }
