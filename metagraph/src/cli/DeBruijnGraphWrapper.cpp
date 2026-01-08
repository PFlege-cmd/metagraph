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
     //TODO: For now,those are hardcoded. get rid of this ugly setup somehow
     unsigned long long num_top_labels = 4294967295;
     const double discovery_fraction = 0.699999999999996;
     const double presence_fraction = 0.0;

     std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>> anchor_coords
                         = this->get_graph()->get_kmer_coordinates(anchor_sequence, num_top_labels,
                                                                   discovery_fraction,
                                                                   presence_fraction);

     std::cout << "Anchor seq: " << anchor_sequence << std::endl;
     //node_index anchor_id = get_graph()->get_graph().find(anchor_sequence);
    //std::cout << "anchor coords Size: " << anchor_coords.size() << std::endl;
     for (unsigned long k = 0; k < anchor_coords.size(); ++k) {
         std::string curr_genome = std::get<0>(anchor_coords[k]);
         //std::cout << "Current genome of coordinate: " << k << ": " << curr_genome << std::endl;
         //std::cout << "Input genome name: " << genome_name << std::endl;
         std::cout << "Curr genome" << curr_genome << std::endl;
         std::cout << "Kmer genome " << std::get<0>(anchor_coords[k]) << std::endl;
         if (std::get<0>(anchor_coords[k]) != genome_name) {
             //throw std::invalid_argument("Anchor not present in genome!");
             continue;
         }


         auto x = std::get<2>(anchor_coords[k])[0];
         //bool found = std::binary_search(x.begin(), x.end(), anchor_position);
         //if (!found)
         //    throw std::invalid_argument("Invalid anchor position");

         // std::cout << "Size of range of coordinates: " <<
         // std::get<2>(edge_node_coords[k])[0].size()<< std::endl;
         auto extracted_coords = std::get<2>(anchor_coords[k]);
         auto coordinate_size = extracted_coords[0].size();
         //std::cout << "Coordinate size is:" << " " << coordinate_size << std::endl;

         for (unsigned long j = 0; j < coordinate_size; ++j) {
             auto all_coords = extracted_coords[0];
             auto potential_anchor_coord = extracted_coords[0][j];
             //std::cout << "Potential coords: " << potential_anchor_coord << std::endl;

             // std::cout << outgoing_edge_coords << std::endl;
             if (potential_anchor_coord == (unsigned long long) anchor_position) {
                 //std::cout << "Found anchor: " << std::endl;
                 /*std::cout << this->get_graph()->get_graph().get_node_sequence(
                         outgoing_nodes[i])
                           << std::endl;*/
                 //std::cout << "Position: " << extracted_coords[0][j] << std::endl;
                 //next_coordinate++; Invalid anchor position
                 // auto next_char = this->get_graph()->get_graph().get_node_sequence(outgoing_nodes[i]).at(size_kmer - 1);
                 /*found_next_coordinate = true;
                 current_kmer = outgoing_nodes[i];
                 outgoing_nodes.clear();*/
                 // starting_kmer.append( 1, next_char);
                 return potential_anchor_coord;
             }
         }
         throw std::invalid_argument("Anchor position not there!");
     }

     throw std::invalid_argument("Invalid anchor position");
}
uint64_t DeBruijnGraphWrapper::retrieveAnchorId(std::string anchor_sequence) {
     std::vector<node_index> nodes;
     nodes.reserve(anchor_sequence.size());
     get_graph()->get_graph().map_to_nodes(anchor_sequence, [&](node_index i) { nodes.push_back(i); });
     //node_index anchor_id = get_graph()->map_to_nodes(get_graph()->get_graph(), anchor_sequence);
     return nodes[0];
}
const char * DeBruijnGraphWrapper::get_sequence_for_coords(std::string genome, unsigned long long start, unsigned long long end, node_index start_index){

//TODO: For now,those are hardcoded. get rid of this ugly setup somehow
     unsigned long long num_top_labels = 4294967295;
     const double discovery_fraction = 0.699999999999996;
     const double presence_fraction = 0.0;

     uint64_t current_kmer = start_index;
     std::vector<uint64_t> outgoing_nodes;
     int size_kmer = this->get_graph()->get_graph().get_k();
     std::string starting_kmer = this->get_graph()->get_graph().get_node_sequence(start_index);

     unsigned long long next_coordinate = start + 1;
     unsigned long long current_coordinate = start;
     //TODO: Update this here correctly , instead of just next_coordinate
     while (current_coordinate != end) {
         bool found_next_coordinate = false;
         bool one_outgoing
                 = this->get_graph()->get_graph().has_single_outgoing(current_kmer);
         bool multi_outgoing
                 = this->get_graph()->get_graph().has_multiple_outgoing(current_kmer);
         std::string one_outgoing_string = (one_outgoing) ? "True" : "False";
         std::string multiple_outgoing_string = (multi_outgoing) ? "True" : "False";
         if (!one_outgoing && !multi_outgoing) {
             std::cout << "Edge stuck!" << std::endl;
             std::cout << next_coordinate << std::endl;
             return "";
         }

         if (multi_outgoing) {
             this->get_graph()->get_graph().adjacent_outgoing_nodes(
                     current_kmer, [&](auto i) { outgoing_nodes.push_back(i); });
             for (unsigned long i = 0; i < outgoing_nodes.size(); ++i) {
                 std::vector<node_index> edge_node = { outgoing_nodes[i] };
                 std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>> edge_node_coords
                         = this->get_graph()->get_kmer_coordinates(edge_node, num_top_labels,
                                                                   discovery_fraction,
                                                                   presence_fraction);
                 for (unsigned long k = 0; k < edge_node_coords.size(); ++k) {
                     std::string curr_genome = std::get<0>(edge_node_coords[k]);
                     if (std::get<0>(edge_node_coords[k]) != genome) {
                         continue;
                     }
                     auto x = std::get<2>(edge_node_coords[k])[0];
                     bool found = std::binary_search(x.begin(), x.end(), next_coordinate);
                     if (!found)
                         continue;

                     auto extracted_coords = std::get<2>(edge_node_coords[k]);
                     auto coordinate_size = extracted_coords[0].size();
                     for (unsigned long j = 0; j < coordinate_size; ++j) {
                         auto all_coords = extracted_coords[0];
                         auto outgoing_edge_coords = extracted_coords[0][j];
                         if (outgoing_edge_coords == next_coordinate) {
                             next_coordinate++;
                             found_next_coordinate = true;
                             current_kmer = outgoing_nodes[i];
                             current_coordinate = outgoing_edge_coords;
                             char next_char = this->get_graph()->get_graph().get_node_sequence(outgoing_nodes[i]).at(size_kmer - 1);
                             starting_kmer.append(1, next_char);
                             outgoing_nodes.clear();
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
             this->get_graph()->get_graph().adjacent_outgoing_nodes(
                     current_kmer, [&](auto i) { outgoing_nodes.push_back(i); });
             for (unsigned long i = 0; i < outgoing_nodes.size(); ++i) {
                 std::vector<node_index> edge_node = { outgoing_nodes[i] };
                 std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>> edge_node_coords
                         = this->get_graph()->get_kmer_coordinates(edge_node, num_top_labels,
                                                                   discovery_fraction,
                                                                   presence_fraction);
                 for (unsigned long k = 0; k < edge_node_coords.size(); ++k) {
                     std::string curr_genome = std::get<0>(edge_node_coords[k]);
                     if (std::get<0>(edge_node_coords[k]) != genome ) {
                         continue;
                     } //TODO: Fix the degenerate issue here!
                     auto x = std::get<2>(edge_node_coords[k])[0];



                     for (unsigned long j = 0;
                          j < std::get<2>(edge_node_coords[k])[0].size(); ++j) {
                         unsigned long long edge_coords = std::get<2>(edge_node_coords[k])[0][j];
                            if (edge_coords == next_coordinate) {
                                found_next_coordinate = true;
                                current_kmer = outgoing_nodes[i];
                                outgoing_nodes.clear();
                                current_coordinate = edge_coords;
                                next_coordinate++;
                                char next_char = this->get_graph()->get_graph().get_node_sequence(outgoing_nodes[i]).at(size_kmer - 1);
                                starting_kmer.append(1, next_char);
                                break;
                            }
                     }
                     if (found_next_coordinate)
                         break;
                 }
                 if (found_next_coordinate)
                     break;
             }
         }
     }

     const  char * result = starting_kmer.c_str();
     return result;
}

uint64_t DeBruijnGraphWrapper::get_first_node_of_coord_range(uint64_t anchor_index, long long start_anchor, std::string genome, long long start_range) {

     //TODO: For now,those are hardcoded. get rid of this ugly setup somehow
     unsigned long long num_top_labels = 4294967295;
     const double discovery_fraction = 0.699999999999996;
     const double presence_fraction = 0.0;

     uint64_t current_kmer = anchor_index;
     std::vector<uint64_t> outgoing_nodes;

     unsigned long long next_coordinate = start_anchor + 1;
     unsigned long long current_coordinate = start_anchor;
     //TODO: Update this here correctly , instead of just next_coordinate
     while (current_coordinate != (unsigned long long) start_range) {
         bool found_next_coordinate = false;
         bool one_outgoing
                 = this->get_graph()->get_graph().has_single_outgoing(current_kmer);
         bool multi_outgoing
                 = this->get_graph()->get_graph().has_multiple_outgoing(current_kmer);
         std::string one_outgoing_string = (one_outgoing) ? "True" : "False";
         std::string multiple_outgoing_string = (multi_outgoing) ? "True" : "False";
         if (!one_outgoing && !multi_outgoing) {
             std::cout << "Edge stuck!" << std::endl;
             std::cout << next_coordinate << std::endl;
             return -1;
         }

         if (multi_outgoing) {
             /*std::cout << "Multiple outgoing in Search for kmer: "
                       << this->get_graph()->get_graph().get_node_sequence(current_kmer)
                       << std::endl;*/
             this->get_graph()->get_graph().adjacent_outgoing_nodes(
                     current_kmer, [&](auto i) { outgoing_nodes.push_back(i); });
             for (unsigned long i = 0; i < outgoing_nodes.size(); ++i) {
                 std::vector<node_index> edge_node = { outgoing_nodes[i] };
                 std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>> edge_node_coords
                         = this->get_graph()->get_kmer_coordinates(edge_node, num_top_labels,
                                                                   discovery_fraction,
                                                                   presence_fraction);
                 for (unsigned long k = 0; k < edge_node_coords.size(); ++k) {
                     std::string curr_genome = std::get<0>(edge_node_coords[k]);
                     if (std::get<0>(edge_node_coords[k]) != genome) {
                         continue;
                     }
                     auto x = std::get<2>(edge_node_coords[k])[0];
                     bool found = std::binary_search(x.begin(), x.end(), next_coordinate);
                     if (!found)
                         continue;

                     auto extracted_coords = std::get<2>(edge_node_coords[k]);
                     auto coordinate_size = extracted_coords[0].size();
                     for (unsigned long j = 0; j < coordinate_size; ++j) {
                         auto all_coords = extracted_coords[0];
                         auto outgoing_edge_coords = extracted_coords[0][j];
                         if (outgoing_edge_coords == next_coordinate) {
                             next_coordinate++;
                             found_next_coordinate = true;
                             current_kmer = outgoing_nodes[i];
                             current_coordinate = outgoing_edge_coords;
                             outgoing_nodes.clear();
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
             this->get_graph()->get_graph().adjacent_outgoing_nodes(
                     current_kmer, [&](auto i) { outgoing_nodes.push_back(i); });
             for (unsigned long i = 0; i < outgoing_nodes.size(); ++i) {
                 std::vector<node_index> edge_node = { outgoing_nodes[i] };
                 std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>> edge_node_coords
                         = this->get_graph()->get_kmer_coordinates(edge_node, num_top_labels,
                                                                   discovery_fraction,
                                                                   presence_fraction);
                 size_t coordinate_number = edge_node_coords.size();

                 for (size_t k = 0; k < coordinate_number; ++k) {
                     std::string curr_genome = std::get<0>(edge_node_coords[k]);
                     if (std::get<0>(edge_node_coords[k]) != genome ) {
                         continue;
                     } //TODO: Fix the degenerate issue here!

                     auto x = std::get<2>(edge_node_coords[k])[0];

                     for (unsigned long j = 0;
                          j < std::get<2>(edge_node_coords[k])[0].size(); ++j) {
                         unsigned long long edge_coords = std::get<2>(edge_node_coords[k])[0][j];

                            if (edge_coords == next_coordinate) {
                                found_next_coordinate = true;
                                current_kmer = outgoing_nodes[i];
                                current_coordinate = edge_coords;
                                outgoing_nodes.clear();
                                next_coordinate++;
                                /*std::cout << "Next coordinate: " << next_coordinate << std::endl;
                                std::cout << "Current coordinate: " << current_coordinate << std::endl;
                                std::cout << "Start: " << start_range << std::endl;*/
                                break;
                            }
                     }
                     if (found_next_coordinate)
                         break;
                 }
                 if (found_next_coordinate)
                     break;
             }
         }
     }
     //std::cout<< "Found first node of range!" << std::endl;
     return current_kmer;
}
