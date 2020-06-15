//
// Created by lessmeaning on 15.06.2020.
//

#include "../../includes/structures/BinaryLifting.h"


int BinaryLifting::get_height(int v, int par) {
    int ret = 0;
    for (int id : tree.g[v]) {
        int to = tree.edges[id].to;
        if (to == par)
            continue;
        ret = std::max(ret, 1 + get_height(to, v));
    }
    return ret;
}

void BinaryLifting::pre_dfs(int v, int par) {
    for (int i = 1; i < LOG; ++i) {
        up[v][i] = up[up[v][i - 1]][i - 1];
        dist[v][i] = dist[v][i - 1] + dist[up[v][i - 1]][i - 1];
    }
    for (int id : tree.g[v]) {
        int to = tree.edges[id].to;
        if (to == par)
            continue;
        up[to][0] = v;
        dist[to][0] = tree.edges[id].weight;
        pre_dfs(to, v);
    }
}

BinaryLifting::BinaryLifting(Tree t) {
    tree = t;
    n = tree.n;
    LOG = log2(get_height(0) + 1) + 2;
    up.resize(n, std::vector<int>(LOG, 0));
    dist.resize(n, std::vector<long long>(LOG, 0));

    pre_dfs(0);
}

BinaryLifting::BinaryLifting() {}