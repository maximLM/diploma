//
// Created by lessmeaning on 15.06.2020.
//

#ifndef DIPLOMA_LADDERDECOMPOSITION_H
#define DIPLOMA_LADDERDECOMPOSITION_H


#include "Tree.h"

class LadderDecomposition {

public:

    Tree tree;
    std::vector<int> height;
    std::vector<int> ind;
    int n;

    void dfs(int v, int par = -1);


    std::vector<int> way;
    std::vector<std::vector<int>> verts;

    LadderDecomposition(Tree t);

};


#endif //DIPLOMA_LADDERDECOMPOSITION_H
