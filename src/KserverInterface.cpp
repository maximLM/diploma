//
// Created by lessmeaning on 15.06.2020.
//


#include "../includes/KserverInterface.h"

KserverInterface::KserverInterface(std::vector<Edge> edgs, std::vector<int> positions) {
    tree = Tree(edgs);
    this->positions = positions;
}

int KserverInterface::serve(int query) {
    assert(0);
}