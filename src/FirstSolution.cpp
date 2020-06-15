//
// Created by lessmeaning on 15.06.2020.
//

#include "../includes/FirstSolution.h"

FirstSolution::FirstSolution(std::vector<Edge> edges, std::vector<int> positions) : KserverInterface(edges, positions) {
    coloringHld = ColoringHLD(Tree(edges));
    movingHld = MovingHLD(Tree(edges));
}

int FirstSolution::serve(int query) {
    int k = positions.size();
    std::vector<int> servers(k);
    for (int i = 0; i < k; ++i) {
        servers[i] = i;
    }
    sort(servers.begin(), servers.end(), [&](int s1, int s2) -> bool {
        int d1 = coloringHld.get_dist(positions[s1], query);
        int d2 = coloringHld.get_dist(positions[s2], query);
        return d1 < d2 || (d1 == d2 && s1 < s2);
    });

    std::vector<int> old_positions = positions;

    coloringHld.color(positions[servers[0]], query, positions[servers[0]] + 1);
    positions[servers[0]] = query;
    for (int i = 1; i < k; ++i) {
        int v = positions[servers[i]];
        pii br = coloringHld.find_color(v, query);
        br.second--;
        int u = movingHld.move(v, query, coloringHld.get_dist(br.first, br.second));
        coloringHld.color(u, v, v + 1);
        positions[servers[i]] = u;
    }

    for (auto v : old_positions) {
        coloringHld.color(v, query, 0);
    }

    return servers[0];
}