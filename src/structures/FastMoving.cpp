//
// Created by lessmeaning on 15.06.2020.
//

#include "../../includes/structures/FastMoving.h"


FastMoving::FastMoving(Tree t): ladderDecomposition(t) {
        tree = Tree(t);
        hld = HLD(tree);
        n = tree.n;
        lifting = BinaryLifting(tree);
}

int FastMoving::move(int from, int to, long long dst) {
    int par = hld.lca(from, to);
    int d = hld.get_dist(from, par);
    if (d < dst) {
        dst -= d;
        dst = hld.get_dist(to, par) - dst;
        from = to;
    }
    int v = from;
    if (!dst)
        return v;
    int j = hld.fastLca.lg2[dst];
    v = lifting.up[v][j];
    dst -= 1 << j;
    return ladderDecomposition.verts[ladderDecomposition.way[v]][ladderDecomposition.ind[v] + dst];
}

long long FastMoving::get_dist(int u, int v) {
    return hld.get_dist(u, v);
}

int FastMoving::lca(int u, int v) {
    return hld.lca(u, v);
}