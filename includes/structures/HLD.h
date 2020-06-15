//
// Created by lessmeaning on 15.06.2020.
//

#ifndef DIPLOMA_HLD_H
#define DIPLOMA_HLD_H

#include "FastLCA.h"

class HLD {
protected:

    Tree tree;
    int n, szway = 0, sztree = 0;
    std::vector<int> way, lf, rf, ind, high, depth, anti_ind;
    std::vector<long long> dist;

    void dfs(int v, int par);

    int pre_dfs(int v, int par);

public:
    FastLCA fastLca;

    HLD(Tree t);

    HLD();

    int lca(int u, int v);

    long long get_dist(int u, int v);
};


#endif //DIPLOMA_HLD_H
