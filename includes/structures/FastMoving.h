//
// Created by lessmeaning on 15.06.2020.
//

#ifndef DIPLOMA_FASTMOVING_H
#define DIPLOMA_FASTMOVING_H


#include "LadderDecomposition.h"
#include "BinaryLifting.h"
#include "HLD.h"
#include "MovingInterface.h"

class FastMoving : public MovingInterface {
    Tree tree;
    BinaryLifting lifting;
    HLD hld;
    LadderDecomposition ladderDecomposition;
    int n;

public:

    FastMoving(Tree t);

    int move(int from, int to, long long dst);

    long long get_dist(int u, int v);

    int lca(int u, int v);
};


#endif //DIPLOMA_FASTMOVING_H
