//
// Created by lessmeaning on 15.06.2020.
//

#include "../includes/ThirdSolution.h"
#include "../includes/structures/FastMoving.h"

ThirdSolution::ThirdSolution(std::vector<Edge> edgs, std::vector<int> positions) : SecondSolution(edgs, positions) {
    FastMoving *tmp = new FastMoving(tree);
    setMovingInstance(tmp);
}