//
// Created by lessmeaning on 15.06.2020.
//

#ifndef DIPLOMA_KSERVERINTERFACE_H
#define DIPLOMA_KSERVERINTERFACE_H

#include <bits/stdc++.h>
#include "structures/Edge.h"
#include "structures/Tree.h"

class KserverInterface {
protected:
    Tree tree;

public:
    std::vector<int> positions;

    KserverInterface(std::vector<Edge> edgs, std::vector<int> positions);

    virtual int serve(int query);
};

#endif //DIPLOMA_KSERVERINTERFACE_H
