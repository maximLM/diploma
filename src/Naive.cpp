//
// Created by lessmeaning on 15.06.2020.
//

#include "../includes/Naive.h"

const int oo = 1e9 + 10;

Naive::Naive(std::vector <Edge> edgs, std::vector<int> positions): KserverInterface(edgs, positions) {}

int Naive::serve(int query) {
    int n = tree.n;
    std::vector<int> dist(n);
    std::vector <std::pair<int, int>> closest(n);
    std::vector<int> server(n, -1);
    for (int i = ((int) positions.size()) - 1; i >= 0; --i) {
        server[positions[i]] = i;
    }

    auto dfs = std::function<void(int, int)>();
    dfs = [&](int v, int par) -> void {
        closest[v] = {oo, oo};
        for (int id : tree.g[v]) {
            int to = tree.edges[id].to;
            if (to == par)
                continue;
            dist[to] = dist[v] + tree.edges[id].weight;
            dfs(to, v);
            closest[v] = min(closest[v], closest[to]);
        }
        if (server[v] != -1)
            closest[v] = {dist[v], server[v]};
    };

    dfs(query, -1);

    auto dfs2 = std::function<void(int, int, long long, int)>();
    dfs2 = [&](int v, int par, long long dst, int cur_server) -> void {
        if (closest[v].first == oo)
            return;
        if (cur_server == closest[v].second) {
            dst = closest[v].first - dist[v];
        } else {
            if (closest[v].first - dist[v] <= dst) {
                cur_server = closest[v].second;
                if (v != query && closest[v].first - dist[par] <= dst)
                    positions[cur_server] = par;
                else
                    positions[cur_server] = v;
                dst = closest[v].first - dist[v];
            }
        }
        for (int id : tree.g[v]) {
            int to = tree.edges[id].to;
            if (to == par)
                continue;
            dfs2(to, v, dst, cur_server);
        }
    };

    dfs2(query, -1, oo, -1);
    return closest[query].second;
}
