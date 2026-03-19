//
// Created by Patrick Flege on 16/03/2026.
//


#include "mock_graph_wrapper.h"

MockGraphWrapper::MockGraphWrapper(mtg::graph::AnnotatedDBG &graph) : GraphWrapper(&graph) {
    //I think that &stuff is converted, ie copied, by copying tyhr pointer to it?

    std::cout << "MockGraphWrapper::MockGraphWrapper()" << &graph  << std::endl;
}

MockGraphWrapper::MockGraphWrapper(MockGraphWrapper& a) : GraphWrapper(a.graph_) {
    this->graph_ = a.graph_;
    std::cout << "MockGraphWrapper::MockGraphWrapper()" << &graph_ << std::endl;
}