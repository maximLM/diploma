//
// Created by lessmeaning on 15.06.2020.
//

#ifndef DIPLOMA_MOVINGINTERFACE_H
#define DIPLOMA_MOVINGINTERFACE_H


class MovingInterface {

public:

    virtual int move(int from, int to, long long dst);

    virtual long long get_dist(int u, int v) ;

    virtual int lca(int u, int v);

    virtual ~MovingInterface();
};


#endif //DIPLOMA_MOVINGINTERFACE_H
