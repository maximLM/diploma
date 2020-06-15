//
// Created by lessmeaning on 15.06.2020.
//

#include "../../includes/structures/LadderDecomposition.h"


void LadderDecomposition::dfs(int v, int par) {
    int mx = -1;
    for (int id : tree.g[v]) {
        int to = tree.edges[id].to;
        if (to == par)
            continue;
        dfs(to, v);
        if (height[to] + 1 > height[v]) {
            height[v] = height[to] + 1;
            mx = to;
        }
    }
    if (mx == -1) {
        way[v] = verts.size();
        verts.push_back({});
    } else {
        way[v] = way[mx];
    }
    ind[v] = verts[way[v]].size();
    verts[way[v]].push_back(v);
}

LadderDecomposition::LadderDecomposition(Tree t) {
    tree = t;
    n = tree.n;
    way.resize(n);
    height.resize(n, 0);
    ind.resize(n);

    dfs(0);
    for (int w = 0; w < verts.size(); ++w) {
        int h = verts[w].size();
        int v = verts[w].back();
        while (h-- || true) {
            if (v == 0)
                break;
            v = tree.parent[v];
            verts[w].push_back(v);
        }
    }
}