//
// Created by lessmeaning on 15.06.2020.
//

#include "Tree.h"

#ifndef DIPLOMA_FASTLCA_H
#define DIPLOMA_FASTLCA_H


class FastLCA {

    Tree tree;
    int n;
    int lg;
    std::vector<int> tin, antin;
    std::vector<std::vector<int>> SP;
    std::vector<int> posl, posr;
    std::vector<int> ord;
    int timer = 0;

    void dfs(int v, int par = -1);

public:

    std::vector<int> lg2;

    FastLCA();

    FastLCA(Tree t);

    int lca(int u, int v);
};


#endif //DIPLOMA_FASTLCA_H
