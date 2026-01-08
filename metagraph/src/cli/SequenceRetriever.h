//
// Created by Patrick Flege on 26/11/2025.
//

#ifndef SEQUENCERETRIEVER_H
#define SEQUENCERETRIEVER_H
#include "GraphWrapper.hpp"


class SequenceRetriever {
    typedef uint64_t node_index;
    public:
        SequenceRetriever(GraphWrapper& wrapper);
        node_index retrieveAnchorId(std::string anchor_sequence);
        node_index retrieveStartId(node_index anchor_index);
        std::string retrieveRegionForRange(node_index start_index);

        void setGenome(const std::string& genome);
        void setAnchorPosition(long long anchor_position);
        void setSequenceStart(long long start);
        void setSequenceEnd(long long end);


    private:
        GraphWrapper& mWrapper;
        std::string mGenome;
        long long mAnchor_position;
        long long mSequenceStart;
        long long mSequenceEnd;
};



#endif //SEQUENCERETRIEVER_H
