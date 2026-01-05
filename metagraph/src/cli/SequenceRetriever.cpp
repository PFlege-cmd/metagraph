//
// Created by Patrick Flege on 26/11/2025.
//

#include "SequenceRetriever.h"

SequenceRetriever::SequenceRetriever(GraphWrapper& wrapper) : mWrapper{wrapper}, mGenome{""}, mAnchor_position{0}, mSequenceStart{0}, mSequenceEnd{0} {
    //std::cout << "Constructor SequenceRetriever" << std::endl;
}

SequenceRetriever::node_index SequenceRetriever::retrieveAnchorId(std::string anchor_sequence) {
    //std::cout << "retrieveanchorId: " << anchor_sequence << std::endl;
    //std::cout << "Address of anchor_sequence: " <<  reinterpret_cast<u_long>(&mWrapper) << std::endl;
    long long retrieved_coord = mWrapper.retrieveAnchorCoordinates(anchor_sequence, mAnchor_position, std::string(mGenome));

    assert(retrieved_coord == this -> mAnchor_position);
    if (retrieved_coord != this -> mAnchor_position) {
        throw std::invalid_argument("Uneqyal pis");
    }
    //std::cout << "retrieved_coord: " << retrieved_coord << std::endl;

    return mWrapper.retrieveAnchorId(anchor_sequence);
}

SequenceRetriever::node_index SequenceRetriever::retrieveStartId(node_index anchor_index) {
    //std::cout << "retrieveStartId: " << anchor_index << std::endl;
    node_index start_node_id = mWrapper.get_first_node_of_coord_range(anchor_index, mAnchor_position, std::string(mGenome), mSequenceStart);
    return start_node_id;
}

std::string SequenceRetriever::retrieveRegionForRange(node_index start_index) {
    std::string retrieved_region = mWrapper.get_sequence_for_coords(mGenome, mSequenceStart, mSequenceEnd, start_index);
    return retrieved_region;
}


void SequenceRetriever::setGenome(std::string genome) {
    mGenome = genome.data();
}

void SequenceRetriever::setSequenceStart(long long sequenceStart) {
    mSequenceStart = sequenceStart;
}

void SequenceRetriever::setSequenceEnd(long long sequenceEnd) {
    mSequenceEnd = sequenceEnd;
}

void SequenceRetriever::setAnchorPosition(long long anchorPosition) {
    mAnchor_position = anchorPosition;
}

