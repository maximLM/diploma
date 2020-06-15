//
// Created by lessmeaning on 15.06.2020.
//

#include "HLD.h"
#include "SegmentTree.h"

#ifndef DIPLOMA_COLORINGHLD_H
#define DIPLOMA_COLORINGHLD_H

typedef std::pair<int, int> pii;

class ColoringHLD : public HLD {
    SegmentTree segmentTree;


    pii find_up(int from, int to);

    pii find_down(int from, int to) ;


public:

    ColoringHLD();

    ColoringHLD(Tree tree);

    void color(int u, int v, int col);

    pii find_color(int from, int to);
};


#endif //DIPLOMA_COLORINGHLD_H
