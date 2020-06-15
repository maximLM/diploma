//
// Created by lessmeaning on 15.06.2020.
//

#include "../../includes/structures/FastLCA.h"

const int oo = 1e9 + 10;

void FastLCA::dfs(int v, int par) {
    tin[v] = timer++;
    antin[tin[v]] = v;
    posr[v] = posl[v] = ord.size();
    ord.push_back(tin[v]);

    for (int id : tree.g[v]) {
        int to = tree.edges[id].to;
        if (to == par)
            continue;
        dfs(to, v);
        posr[v] = ord.size();
        ord.push_back(tin[v]);
    }
}

FastLCA::FastLCA() {}

FastLCA::FastLCA(Tree t) {
    tree = t;
    n = tree.n;
    tin.resize(n);
    antin.resize(n);
    lg = log2(n) + 2;
    posl.resize(n);
    posr.resize(n);

    dfs(0);

    SP.resize(lg, std::vector<int>(ord.size(), oo));

    SP[0] = ord;

    for (int j = 1; j < lg; ++j) {
        for (int i = 0; i < ord.size(); ++i) {
            SP[j][i] = SP[j - 1][i];
            int l = i + (1 << j - 1);
            if (l < ord.size())
                SP[j][i] = std::min(SP[j][i], SP[j - 1][l]);
        }
    }
    lg2.resize(ord.size() + 1);
    int j = 0;
    lg2[0] = j;
    for (int i = 1; i < ord.size() + 1; ++i) {
        while ((1 << j) <= i)
            ++j;
        --j;
        lg2[i] = j;
    }
}

int FastLCA::lca(int u, int v) {
    if (tin[u] > tin[v])
        std::swap(u, v);
    int l = posl[u];
    int r = posr[v];
    int len = r - l + 1;
    int j = lg2[len];
    int tim = std::min(SP[j][l], SP[j][r - (1 << j) + 1]);
    return antin[tim];
}