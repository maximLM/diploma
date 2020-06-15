//
// Created by lessmeaning on 15.06.2020.
//

#ifndef DIPLOMA_NAIVE_H
#define DIPLOMA_NAIVE_H

#include "KserverInterface.h"
#include "structures/Tree.h"

class Naive : public KserverInterface {

public:


    Naive(std::vector<Edge> edgs, std::vector<int> positions);

    int serve(int query);
};

#endif //DIPLOMA_NAIVE_H
