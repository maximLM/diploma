//
// Created by lessmeaning on 15.06.2020.
//

#ifndef DIPLOMA_BINARYLIFTING_H
#define DIPLOMA_BINARYLIFTING_H

#include "Tree.h"

class BinaryLifting {

    Tree tree;
    std::vector<std::vector<long long>> dist;
    int n, LOG;

    int get_height(int v, int par = -1);

    void pre_dfs(int v, int par = -1);

public:

    std::vector<std::vector<int>> up;

    BinaryLifting(Tree t);

    BinaryLifting();
};


#endif //DIPLOMA_BINARYLIFTING_H
