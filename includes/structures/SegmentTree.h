//
// Created by lessmeaning on 15.06.2020.
//
#include <bits/stdc++.h>

#ifndef DIPLOMA_SEGMENTTREE_H
#define DIPLOMA_SEGMENTTREE_H

typedef std::pair<int, int> pii;

class SegmentTree {
    std::vector<int> t;
    std::vector<int> pushing;
    int n;

    void push(int v);

    void upd(int v, int tl, int tr, int l, int r, int col);

    pii get_left(int v, int tl, int tr, int l, int r);

    pii get_right(int v, int tl, int tr, int l, int r);

public:

    SegmentTree(int n) ;

    SegmentTree() ;

    void upd(int l, int r, int col);

    pii get_left(int l, int r);

    pii get_right(int l, int r) ;

};


#endif //DIPLOMA_SEGMENTTREE_H
