//
// Created by lessmeaning on 15.06.2020.
//

#include "../../includes/structures/Tree.h"


void Tree::dfs(int v, int par) {
    parent[v] = par;
    for (int id : g[v]) {
        int to = edges[id].to;
        if (to == par)
            continue;
        parent_edge[to] = id;
        dfs(to, v);
    }
}

Tree::Tree(std::vector <Edge> eds, int rt) {
    root = rt;
    n = eds.size() + 1;
    g.resize(n);
    parent.resize(n);
    parent_edge.resize(n, -1);
    for (auto e : eds) {
        g[e.from].push_back(edges.size());
        edges.push_back(e);
        std::swap(e.from, e.to);
        g[e.from].push_back(edges.size());
        edges.push_back(e);
    }
    dfs(root);
}

Tree::Tree() {}