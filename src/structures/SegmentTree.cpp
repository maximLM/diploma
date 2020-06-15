//
// Created by lessmeaning on 15.06.2020.
//

#include "../../includes/structures/SegmentTree.h"

void SegmentTree::push(int v) {
    if (pushing[v] != -1) {
        pushing[v << 1] = pushing[v];
        pushing[v << 1 | 1] = pushing[v];
        t[v << 1] = pushing[v];
        t[v << 1 | 1] = pushing[v];
        pushing[v] = -1;
    }
}

void SegmentTree::upd(int v, int tl, int tr, int l, int r, int col) {
    if (l > r)
        return;
    if (tl == l && tr == r) {
        t[v] = col;
        pushing[v] = col;
        return;
    }
    push(v);
    int tm = tl + tr >> 1;
    upd(v << 1, tl, tm, l, std::min(tm, r), col);
    upd(v << 1 | 1, tm + 1, tr, std::max(tm + 1, l), r, col);
    t[v] = std::max(t[v << 1], t[v << 1 | 1]);
}

pii SegmentTree::get_left(int v, int tl, int tr, int l, int r) {
    if (l > r || !t[v])
        return {0, 0};
    if (tl == tr)
        return {tl, t[v]};
    push(v);
    int tm = tl + tr >> 1;
    pii val = get_left(v << 1, tl, tm, l, std::min(tm, r));
    if (!val.second)
        val = get_left(v << 1 | 1, tm + 1, tr, std::max(tm + 1, l), r);
    return val;
}

pii SegmentTree::get_right(int v, int tl, int tr, int l, int r) {
    if (l > r || !t[v])
        return {0, 0};
    if (tl == tr)
        return {tl, t[v]};
    push(v);
    int tm = tl + tr >> 1;
    pii val = get_right(v << 1 | 1, tm + 1, tr, std::max(tm + 1, l), r);
    if (!val.second)
        val = get_right(v << 1, tl, tm, l, std::min(tm, r));
    return val;
}


SegmentTree::SegmentTree(int n) {
    this->n = n;
    t.resize(n << 2, 0);
    pushing.resize(n << 2, -1);
}

SegmentTree::SegmentTree() {
    // dummy
}

void SegmentTree::upd(int l, int r, int col) {
    upd(1, 0, n - 1, l, r, col);
}

pii SegmentTree::get_left(int l, int r) {
    return get_left(1, 0, n - 1, l, r);
}

pii SegmentTree::get_right(int l, int r) {
    return get_right(1, 0, n - 1, l, r);
}
