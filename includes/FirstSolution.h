//
// Created by lessmeaning on 15.06.2020.
//

#ifndef DIPLOMA_FIRSTSOLUTION_H
#define DIPLOMA_FIRSTSOLUTION_H


#include "KserverInterface.h"
#include "structures/ColoringHLD.h"
#include "structures/MovingHLD.h"

class FirstSolution : public KserverInterface {
    ColoringHLD coloringHld;
    MovingHLD movingHld;

public:

    FirstSolution(std::vector<Edge> edges, std::vector<int> positions);

    int serve(int query);
};


#endif //DIPLOMA_FIRSTSOLUTION_H
