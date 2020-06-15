//
// Created by lessmeaning on 15.06.2020.
//

#include "../../includes/structures/MovingHLD.h"

MovingHLD::MovingHLD() : HLD() {}

MovingHLD::MovingHLD(Tree tree) : HLD(tree) {}

int MovingHLD::move_up(int from, int to, long long &dst) {
    while (way[from] != way[to]) {
        int wind = way[from];
        if (dst >= dist[from] - dist[high[wind]]) {
            dst -= dist[from] - dist[high[wind]];
        } else {
            break;
        }
        from = high[wind];
        if (dst >= tree.edges[tree.parent_edge[from]].weight) {
            dst -= tree.edges[tree.parent_edge[from]].weight;
        } else {
            break;
        }
        from = tree.parent[from];
    }
    int wind = way[from];
    int l = lf[wind] + ind[from];
    int r = wind == way[to] ? lf[wind] + ind[to] + 1 : rf[wind] + 1;
    while (r - l > 1) {
        int m = r + l >> 1;
        if (dist[from] - dist[anti_ind[m]] <= dst) {
            l = m;
        } else {
            r = m;
        }
    }
    int ret = anti_ind[l];
    dst -= dist[from] - dist[ret];
    return ret;
}

int MovingHLD::move_down(int from, int to, long long &dst) {
    std::vector<int> highs;
    int v = to;
    while (way[from] != way[v]) {
        highs.push_back(high[way[v]]);
        v = tree.parent[high[way[v]]];
    }

    v = from;
    while (way[v] != way[to]) {
        int nxt = highs.back();

        int u = tree.parent[nxt];
        if (dst >= dist[u] - dist[v]) {
            dst -= dist[u] - dist[v];
        } else {
            break;
        }
        v = u;

        if (dst >= tree.edges[tree.parent_edge[nxt]].weight) {
            dst -= tree.edges[tree.parent_edge[nxt]].weight;
        } else {
            break;
        }
        v = nxt;
        highs.pop_back();
    }
    int r = lf[way[v]] + ind[v];
    int l = way[v] == way[to] ? lf[way[to]] + ind[to] - 1 : lf[way[v]] + ind[tree.parent[highs.back()]] - 1;
    while (r - l > 1) {
        int m = r + l >> 1;
        if (dist[anti_ind[m]] - dist[v] <= dst) {
            r = m;
        } else {
            l = m;
        }
    }
    int ret = anti_ind[r];
    dst -= dist[ret] - dist[v];
    return ret;
}

int MovingHLD::move(int from, int to, long long dst) {
    int par = HLD::lca(from, to);
    int ans = move_up(from, par, dst);
    if (ans == par)
        ans = move_down(par, to, dst);
    return ans;
}

int MovingHLD::lca(int u, int v) {
    return HLD::lca(u, v);
}

long long MovingHLD::get_dist(int u, int v) {
    return HLD::get_dist(u, v);
}