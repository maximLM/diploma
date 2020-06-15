//
// Created by lessmeaning on 15.06.2020.
//

#ifndef DIPLOMA_MOVINGHLD_H
#define DIPLOMA_MOVINGHLD_H


#include "HLD.h"
#include "MovingInterface.h"

class MovingHLD : HLD, public MovingInterface{

    int move_up(int from, int to, long long & dst);

    int move_down(int from, int to, long long & dst);

public:

    MovingHLD();

    MovingHLD(Tree tree);

    int move(int from, int to, long long dst);

    int lca(int u, int v) ;

    long long get_dist(int u, int v);
};


#endif //DIPLOMA_MOVINGHLD_H
