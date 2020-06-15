//
// Created by lessmeaning on 15.06.2020.
//
#include<bits/stdc++.h>
#include "Edge.h"

#ifndef DIPLOMA_TREE_H
#define DIPLOMA_TREE_H


class Tree {

    void dfs(int v, int par = -1);

public:

    int n;
    std::vector<Edge> edges;
    std::vector<std::vector<int>> g;
    std::vector<int> parent, parent_edge;
    int root;

    Tree(std::vector<Edge> eds, int rt = 0);

    Tree();

    int get_size();
};


#endif //DIPLOMA_TREE_H
