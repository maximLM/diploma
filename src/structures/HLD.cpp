//
// Created by lessmeaning on 15.06.2020.
//

#include "../../includes/structures/HLD.h"

void HLD::dfs(int v, int par) {
    bool leaf = true;
    for (auto id : tree.g[v]) {
        int to = tree.edges[id].to;
        if (to == par)
            continue;
        leaf = false;
    }
    if (leaf) {
        way[v] = szway++;
        int wind = way[v];
        rf[wind] = lf[wind] = sztree++;
        ind[v] = 0;
        high[wind] = v;
        anti_ind[ind[v] + lf[wind]] = v;
        return;
    }
    for (int i = 1; i < (int) tree.g[v].size(); ++i) {
        int id = tree.g[v][i];
        int to = tree.edges[id].to;
        if (to == par)
            continue;
        dfs(to, v);
    }

    int to = tree.edges[tree.g[v][0]].to;
    dfs(to, v);
    int wind = way[to];
    way[v] = wind;
    rf[wind] = sztree++;
    ind[v] = ind[to] + 1;
    anti_ind[ind[v] + lf[wind]] = v;
    high[wind] = v;
}

int HLD::pre_dfs(int v, int par) {
    depth[v] = par == -1 ? 0 : depth[par] + 1;
    int cnt = 1;
    int mx = -1;
    int idmx = -1;
    for (int i = 0; i < tree.g[v].size(); ++i) {
        int id = tree.g[v][i];
        int to = tree.edges[id].to;
        if (to == par)
            continue;

        dist[to] = dist[v] + tree.edges[id].weight;
        int cur = pre_dfs(to, v);
        if (cur > mx) {
            mx = cur;
            idmx = i;
        }
        cnt += cur;
    }
    if (idmx != -1) {
        std::swap(tree.g[v][0], tree.g[v][idmx]);
    }
    return cnt;
}

HLD::HLD(Tree t) {
    tree = t;
    fastLca = FastLCA(tree);
    n = t.n;
    way.resize(n, 0);
    lf.resize(n, 0);
    rf.resize(n, 0);
    ind.resize(n, 0);
    high.resize(n, 0);
    depth.resize(n, 0);
    dist.resize(n, 0);
    anti_ind.resize(n, 0);

    pre_dfs(0, -1);
    dfs(0, -1);
}

HLD::HLD() {}

int HLD::lca(int u, int v) {
    return fastLca.lca(u, v);
}

long long HLD::get_dist(int u, int v) {
    return dist[u] + dist[v] - dist[lca(u, v)] * 2;
}