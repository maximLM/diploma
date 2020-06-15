//
// Created by lessmeaning on 15.06.2020.
//

#ifndef DIPLOMA_SECONDSOLUTION_H
#define DIPLOMA_SECONDSOLUTION_H


#include "KserverInterface.h"
#include "structures/MovingInterface.h"

class SecondSolution: public KserverInterface {
protected:

    MovingInterface * movingInstance;
    Tree tree;
    std::vector<int> tin, tout, server;
    std::vector<std::vector<int>> g;
    std::vector<long long> dist;
    std::vector<std::pair<long long, int>> closest;
    int n;

    int timer = 1;

    void dfs(int v, int par = -1);

public:

    SecondSolution(std::vector<Edge> edges, std::vector<int> positions);

    void setMovingInstance(MovingInterface * moving);

    ~SecondSolution();


    bool is_parent(int u, int v);

    int serve(int query);
};


#endif //DIPLOMA_SECONDSOLUTION_H
