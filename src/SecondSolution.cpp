//
// Created by lessmeaning on 15.06.2020.
//

#include "../includes/SecondSolution.h"
#include "../includes/structures/MovingHLD.h"

const int oo = 1e9 + 10;

void SecondSolution::dfs(int v, int par) {
    tin[v] = timer++;
    tout[v] = tin[v];
    for (int id : tree.g[v]) {
        int to = tree.edges[id].to;
        if (to == par)
            continue;
        dfs(to, v);
        tout[v] = tout[to];
    }
}

SecondSolution::SecondSolution(std::vector<Edge> edges, std::vector<int> positions): KserverInterface(edges, positions) {
    tree = Tree(edges);
    this->positions = positions;
    n = tree.n;

    MovingHLD *tmp = new MovingHLD(tree);
    movingInstance = tmp;

    dist.resize(n);
    tin.resize(n);
    tout.resize(n);
    g.resize(n);
    closest.resize(n);
    server.resize(n, -1);

    dfs(0);
}

void SecondSolution::setMovingInstance(MovingInterface *moving) {
    delete movingInstance;
    movingInstance = moving;
}

bool SecondSolution::is_parent(int u, int v) {
    return tin[u] <= tin[v] && tin[v] <= tout[u];
}

int SecondSolution::serve(int query) {
    for (int i = ((int) positions.size()) - 1; i >= 0; --i) {
        server[positions[i]] = i;
    }
    std::vector<int> verts = positions;
    verts.push_back(query);
    auto cmp_less = [&](int v1, int v2) -> bool {
        return tin[v1] < tin[v2];
    };
    auto cmp_equal = [&](int v1, int v2) -> bool {
        return tin[v1] == tin[v2];
    };
    sort(verts.begin(), verts.end(), cmp_less);
    verts.resize(unique(verts.begin(), verts.end(), cmp_equal) - verts.begin());

    std::vector<int> tmp = verts;
    for (int i = 1; i < verts.size(); ++i) {
        tmp.push_back(movingInstance->lca(verts[i - 1], verts[i]));
    }

    verts = tmp;
    sort(verts.begin(), verts.end(), cmp_less);
    verts.resize(unique(verts.begin(), verts.end(), cmp_equal) - verts.begin());

    std::vector<int> stack;
    std::vector<Edge> edges;

    for (auto v : verts) {
        while (!stack.empty() && !is_parent(stack.back(), v)) {
            stack.pop_back();
        }
        if (!stack.empty()) {
            g[stack.back()].push_back(edges.size());
            long long len = movingInstance->get_dist(stack.back(), v);
            edges.push_back({stack.back(), v, len});
            g[v].push_back(edges.size());
            edges.push_back({v, stack.back(), len});
        }
        stack.push_back(v);
    }

    auto dfs = std::function<void(int, int)>();
    dfs = [&](int v, int par) -> void {
        closest[v] = {oo, oo};
        for (int id : g[v]) {
            int to = edges[id].to;
            if (to == par)
                continue;
            dist[to] = dist[v] + edges[id].weight;
            dfs(to, v);
            closest[v] = min(closest[v], closest[to]);
        }
        if (server[v] != -1) {
            closest[v] = {dist[v], server[v]};
        }
    };
    dist[query] = 0;
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
                if (v != query) {
                    dst -= closest[v].first - dist[v];
                    positions[cur_server] = movingInstance->move(v, query, dst);
                } else {
                    positions[cur_server] = v;
                }
                dst = closest[v].first - dist[v];
            }
        }
        for (int id : g[v]) {
            int to = edges[id].to;
            if (to == par)
                continue;
            dfs2(to, v, dst, cur_server);
        }
    };

    dfs2(query, -1, oo, -1);
    int ret = closest[query].second;
    for (auto v : verts) {
        server[v] = -1;
        g[v].clear();
    }

    return ret;
}

SecondSolution::~SecondSolution() {
    delete movingInstance;
}
