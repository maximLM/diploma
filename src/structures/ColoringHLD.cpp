//
// Created by lessmeaning on 15.06.2020.
//

#include "../../includes/structures/ColoringHLD.h"


ColoringHLD::ColoringHLD() : HLD() {}

ColoringHLD::ColoringHLD(Tree tree) : HLD(tree) {
        segmentTree = SegmentTree(n);
}

void ColoringHLD::color(int u, int v, int col) {
    while (way[u] != way[v]) {
        if (depth[high[way[u]]] > depth[high[way[v]]]) // u is now higher
            std::swap(u, v);
        int wind = way[v];
        segmentTree.upd(lf[wind] + ind[v], rf[wind], col);
        v = tree.parent[high[wind]];
    }
    if (depth[u] > depth[v]) // u is now higher
        std::swap(u, v);
    int wind = way[v];
    segmentTree.upd(lf[wind] + ind[v], lf[wind] + ind[u], col);
}

pii ColoringHLD::find_up(int from, int to) {
    while (way[from] != way[to]) {
        int wind = way[from];
        pii cur = segmentTree.get_left(lf[wind] + ind[from], rf[wind]);
        if (cur.second)
            return cur;
        from = tree.parent[high[wind]];
    }
    int wind = way[from];
    return segmentTree.get_left(lf[wind] + ind[from], lf[wind] + ind[to]);
}

pii ColoringHLD::find_down(int from, int to) {
    pii ret = {0, 0};
    while (way[to] != way[from]) {
        int wind = way[to];
        pii cur = segmentTree.get_right(lf[wind] + ind[to], rf[wind]);
        if (cur.second)
            ret = cur;
        to = tree.parent[high[wind]];
    }

    int wind = way[to];
    pii cur = segmentTree.get_right(lf[wind] + ind[to], lf[wind] + ind[from]);
    if (cur.second)
        ret = cur;
    return ret;
}

pii ColoringHLD::find_color(int from, int to) {
    int par = lca(from, to);
    pii ret = find_up(from, par);
    if (!ret.second)
        ret = find_down(par, to);
    assert(ret.second);
    ret.first = anti_ind[ret.first];
    return ret;
}