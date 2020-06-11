#include <bits/stdc++.h>

#define watch(x) cout << (#x) << " = " << (x) << endl;

using namespace std;


typedef long long ll;
typedef pair<int, int> pii;

const ll oo = 1e9 + 10;

struct Edge {
    int from, to;
    ll weight;
};

struct Tree {
    int n;
    vector<Edge> edges;
    vector<vector<int>> g;
    vector<int> parent, parent_edge;
    int root;

    void dfs(int v, int par = -1) {
        parent[v] = par;
        for (int id : g[v]) {
            int to = edges[id].to;
            if (to == par)
                continue;
            parent_edge[to] = id;
            dfs(to, v);
        }
    }

    Tree(vector<Edge> eds, int rt = 0) {
        root = rt;
        n = eds.size() + 1;
        g.resize(n);
        parent.resize(n);
        parent_edge.resize(n, -1);
        for (auto e : eds) {
            g[e.from].push_back(edges.size());
            edges.push_back(e);
            swap(e.from, e.to);
            g[e.from].push_back(edges.size());
            edges.push_back(e);
        }
        dfs(root);
    }

    Tree() {}

};

struct SegmentTree {
    vector<int> t;
    vector<int> pushing;
    int n;


private:

    void push(int v) {
        if (pushing[v] != -1) {
            pushing[v << 1] = pushing[v];
            pushing[v << 1 | 1] = pushing[v];
            t[v << 1] = pushing[v];
            t[v << 1 | 1] = pushing[v];
            pushing[v] = -1;
        }
    }

    void upd(int v, int tl, int tr, int l, int r, int col) {
        if (l > r)
            return;
        if (tl == l && tr == r) {
            t[v] = col;
            pushing[v] = col;
            return;
        }
        push(v);
        int tm = tl + tr >> 1;
        upd(v << 1, tl, tm, l, min(tm, r), col);
        upd(v << 1 | 1, tm + 1, tr, max(tm + 1, l), r, col);
        t[v] = max(t[v << 1], t[v << 1 | 1]);
    }

    pii get_left(int v, int tl, int tr, int l, int r) {
        if (l > r || !t[v])
            return {0, 0};
        if (tl == tr)
            return {tl, t[v]};
        push(v);
        int tm = tl + tr >> 1;
        pii val = get_left(v << 1, tl, tm, l, min(tm, r));
        if (!val.second)
            val = get_left(v << 1 | 1, tm + 1, tr, max(tm + 1, l), r);
        return val;
    }

    pii get_right(int v, int tl, int tr, int l, int r) {
        if (l > r || !t[v])
            return {0, 0};
        if (tl == tr)
            return {tl, t[v]};
        push(v);
        int tm = tl + tr >> 1;
        pii val = get_right(v << 1 | 1, tm + 1, tr, max(tm + 1, l), r);
        if (!val.second)
            val = get_right(v << 1, tl, tm, l, min(tm, r));
        return val;
    }

public:

    SegmentTree(int n) {
        this->n = n;
        t.resize(n << 2, 0);
        pushing.resize(n << 2, -1);
    }

    SegmentTree() {
        // dummy
    }

    void upd(int l, int r, int col) {
        upd(1, 0, n - 1, l, r, col);
    }

    pii get_left(int l, int r) {
        return get_left(1, 0, n - 1, l, r);
    }

    pii get_right(int l, int r) {
        return get_right(1, 0, n - 1, l, r);
    }

};

struct HLD {
    Tree tree;
    int n, szway = 0, sztree = 0;
    vector<int> way, lf, rf, ind, high, depth, anti_ind;
    vector<ll> dist;

    void dfs(int v, int par) {
        bool leaf = true;
        for (auto id : tree.g[v]) {
            int to = tree.edges[id].to;
            if (to == par)
                continue;
            leaf = false;
        }
        if (leaf) {
            way[v] = szway++;
            int wind = way[v];
            rf[wind] = lf[wind] = sztree++;
            ind[v] = 0;
            high[wind] = v;
            anti_ind[ind[v] + lf[wind]] = v;
            return;
        }
        for (int i = 1; i < (int) tree.g[v].size(); ++i) {
            int id = tree.g[v][i];
            int to = tree.edges[id].to;
            if (to == par)
                continue;
            dfs(to, v);
        }

        int to = tree.edges[tree.g[v][0]].to;
        dfs(to, v);
        int wind = way[to];
        way[v] = wind;
        rf[wind] = sztree++;
        ind[v] = ind[to] + 1;
        anti_ind[ind[v] + lf[wind]] = v;
        high[wind] = v;
    }

    int pre_dfs(int v, int par) {
        depth[v] = par == -1 ? 0 : depth[par] + 1;
        int cnt = 1;
        int mx = -1;
        int idmx = -1;
        for (int i = 0; i < tree.g[v].size(); ++i) {
            int id = tree.g[v][i];
            int to = tree.edges[id].to;
            if (to == par)
                continue;

            dist[to] = dist[v] + tree.edges[id].weight;
            int cur = pre_dfs(to, v);
            if (cur > mx) {
                mx = cur;
                idmx = i;
            }
            cnt += cur;
        }
        if (idmx != -1) {
            swap(tree.g[v][0], tree.g[v][idmx]);
        }
        return cnt;
    }

    HLD(Tree t) {
        tree = t;
        n = t.n;
        way.resize(n, 0);
        lf.resize(n, 0);
        rf.resize(n, 0);
        ind.resize(n, 0);
        high.resize(n, 0);
        depth.resize(n, 0);
        dist.resize(n, 0);
        anti_ind.resize(n, 0);

        pre_dfs(0, -1);
        dfs(0, -1);
    }

    HLD() {}

    int lca(int u, int v) {
        while (way[u] != way[v]) {
            if (depth[high[way[u]]] > depth[high[way[v]]]) // u is now higher
                swap(u, v);
            v = tree.parent[high[way[v]]];
        }
        return depth[u] > depth[v] ? v : u;
    }

    ll get_dist(int u, int v) {
        return dist[u] + dist[v] - dist[lca(u, v)] * 2;
    }

};

struct ColoringHLD : HLD {

    SegmentTree segmentTree;

    ColoringHLD() : HLD() {}

    ColoringHLD(Tree tree) : HLD(tree) {
        segmentTree = SegmentTree(n);
    }

    void color(int u, int v, int col) {
        while (way[u] != way[v]) {
            if (depth[high[way[u]]] > depth[high[way[v]]]) // u is now higher
                swap(u, v);
            int wind = way[v];
            segmentTree.upd(lf[wind] + ind[v], rf[wind], col);
            v = tree.parent[high[wind]];
        }
        if (depth[u] > depth[v]) // u is now higher
            swap(u, v);
        int wind = way[v];
        segmentTree.upd(lf[wind] + ind[v], lf[wind] + ind[u], col);
    }

    pii find_up(int from, int to) {
        while (way[from] != way[to]) {
            int wind = way[from];
            pii cur = segmentTree.get_left(lf[wind] + ind[from], rf[wind]);
            if (cur.second)
                return cur;
            from = tree.parent[high[wind]];
        }
        int wind = way[from];
        return segmentTree.get_left(lf[wind] + ind[from], lf[wind] + ind[to]);
    }

    pii find_down(int from, int to) {
        pii ret = {0, 0};
        while (way[to] != way[from]) {
            int wind = way[to];
            pii cur = segmentTree.get_right(lf[wind] + ind[to], rf[wind]);
            if (cur.second)
                ret = cur;
            to = tree.parent[high[wind]];
        }

        int wind = way[to];
        pii cur = segmentTree.get_right(lf[wind] + ind[to], lf[wind] + ind[from]);
        if (cur.second)
            ret = cur;
        return ret;
    }

    pii find_color(int from, int to) {
        int par = lca(from, to);
        pii ret = find_up(from, par);
        if (!ret.second)
            ret = find_down(par, to);
        assert(ret.second);
        ret.first = anti_ind[ret.first];
        return ret;
    }
};

struct MovingInterface {
    virtual int move(int from, int to, ll dst) {
        assert(0);
        return -1;
    }
    virtual ll get_dist(int u, int v) {
        assert(0);
        return -1;
    }
    virtual int lca(int u, int v) {
        assert(0);
        return -1;
    }

    virtual ~MovingInterface() {}
};

struct MovingHLD : HLD, MovingInterface {


    MovingHLD() : HLD() {}

    MovingHLD(Tree tree) : HLD(tree) {}

    int move_up(int from, int to, ll &dst) {
        while (way[from] != way[to]) {
            int wind = way[from];
            if (dst >= dist[from] - dist[high[wind]]) {
                dst -= dist[from] - dist[high[wind]];
            } else {
                break;
            }
            from = high[wind];
            if (dst >= tree.edges[tree.parent_edge[from]].weight) {
                dst -= tree.edges[tree.parent_edge[from]].weight;
            } else {
                break;
            }
            from = tree.parent[from];
        }
        int wind = way[from];
        int l = lf[wind] + ind[from];
        int r = wind == way[to] ? lf[wind] + ind[to] + 1 : rf[wind] + 1;
        while (r - l > 1) {
            int m = r + l >> 1;
            if (dist[from] - dist[anti_ind[m]] <= dst) {
                l = m;
            } else {
                r = m;
            }
        }
        int ret = anti_ind[l];
        dst -= dist[from] - dist[ret];
        return ret;
    }

    int move_down(int from, int to, ll &dst) {
        vector<int> highs;
        int v = to;
        while (way[from] != way[v]) {
            highs.push_back(high[way[v]]);
            v = tree.parent[high[way[v]]];
        }

        v = from;
        while (way[v] != way[to]) {
            int nxt = highs.back();

            int u = tree.parent[nxt];
            if (dst >= dist[u] - dist[v]) {
                dst -= dist[u] - dist[v];
            } else {
                break;
            }
            v = u;

            if (dst >= tree.edges[tree.parent_edge[nxt]].weight) {
                dst -= tree.edges[tree.parent_edge[nxt]].weight;
            } else {
                break;
            }
            v = nxt;
            highs.pop_back();
        }
        int r = lf[way[v]] + ind[v];
        int l = way[v] == way[to] ? lf[way[to]] + ind[to] - 1 : lf[way[v]] + ind[tree.parent[highs.back()]] - 1;
        while (r - l > 1) {
            int m = r + l >> 1;
            if (dist[anti_ind[m]] - dist[v] <= dst) {
                r = m;
            } else {
                l = m;
            }
        }
        int ret = anti_ind[r];
        dst -= dist[ret] - dist[v];
        return ret;
    }

    int move(int from, int to, ll dst) {
        int par = HLD::lca(from, to);
        int ans = move_up(from, par, dst);
        if (ans == par)
            ans = move_down(par, to, dst);
        return ans;
    }

    int lca(int u, int v) {
        return HLD::lca(u, v);
    }

    ll get_dist(int u, int v) {
        return HLD::get_dist(u, v);
    }

};

struct KServers {
    ColoringHLD coloringHld;
    MovingHLD movingHld;
    vector<int> positions;

    KServers(vector<Edge> edges, vector<int> positions) {
        coloringHld = ColoringHLD(Tree(edges));
        movingHld = MovingHLD(Tree(edges));
        this->positions = positions;
    }

    int serve(int query) {
        int k = positions.size();
        vector<int> servers(k);
        for (int i = 0; i < k; ++i) {
            servers[i] = i;
        }
        sort(servers.begin(), servers.end(), [&](int s1, int s2) -> bool {
            int d1 = coloringHld.get_dist(positions[s1], query);
            int d2 = coloringHld.get_dist(positions[s2], query);
            return d1 < d2 || (d1 == d2 && s1 < s2);
        });

        vector<int> old_positions = positions;

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
};

struct Naive {
    Tree tree;
    vector<int> positions;

    Naive(vector<Edge> edgs, vector<int> positions) {
        this->positions = positions;
        tree = Tree(edgs);
    }

    int serve(int query) {
        int n = tree.n;
        vector<int> dist(n);
        vector<pii> closest(n);
        vector<int> server(n, -1);
        for (int i = ((int) positions.size()) - 1; i >= 0; --i) {
            server[positions[i]] = i;
        }

        auto dfs = function<void(int, int)>();
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

        auto dfs2 = function<void(int, int, ll, int)>();
        dfs2 = [&](int v, int par, ll dst, int cur_server) -> void {
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
};

struct FasterKServers {
    MovingInterface * movingInstance;
    Tree tree;
    vector<int> positions, tin, tout, server;
    vector<vector<int>> g;
    vector<ll> dist;
    vector<pair<ll, int>> closest;
    int n;

    int timer = 1;

    void dfs(int v, int par = -1) {
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

    FasterKServers(vector<Edge> edges, vector<int> positions) {
        tree = Tree(edges);
        this->positions = positions;
        n = tree.n;

        MovingHLD * tmp = new MovingHLD(tree);
        movingInstance = tmp;

        dist.resize(n);
        tin.resize(n);
        tout.resize(n);
        g.resize(n);
        closest.resize(n);
        server.resize(n, -1);

        dfs(0);
    }

    void setMovingInstance(MovingInterface * moving) {
        delete movingInstance;
        movingInstance = moving;
    }

    bool is_parent(int u, int v) {
        return tin[u] <= tin[v] && tin[v] <= tout[u];
    }

    int serve(int query) {
        for (int i = ((int) positions.size()) - 1; i >= 0; --i) {
            server[positions[i]] = i;
        }
        vector<int> verts = positions;
        verts.push_back(query);
        auto cmp_less = [&](int v1, int v2) -> bool {
            return tin[v1] < tin[v2];
        };
        auto cmp_equal = [&](int v1, int v2) -> bool {
            return tin[v1] == tin[v2];
        };
        sort(verts.begin(), verts.end(), cmp_less);
        verts.resize(unique(verts.begin(), verts.end(), cmp_equal) - verts.begin());

        vector<int> tmp = verts;
        for (int i = 1; i < verts.size(); ++i) {
            tmp.push_back(movingInstance->lca(verts[i - 1], verts[i]));
        }

        verts = tmp;
        sort(verts.begin(), verts.end(), cmp_less);
        verts.resize(unique(verts.begin(), verts.end(), cmp_equal) - verts.begin());

        vector<int> stack;
        vector<Edge> edges;

        for (auto v : verts) {
            while (!stack.empty() && !is_parent(stack.back(), v)) {
                stack.pop_back();
            }
            if (!stack.empty()) {
                g[stack.back()].push_back(edges.size());
                ll len = movingInstance->get_dist(stack.back(), v);
                edges.push_back({stack.back(), v, len});
                g[v].push_back(edges.size());
                edges.push_back({v, stack.back(), len});
            }
            stack.push_back(v);
        }

        auto dfs = function<void(int, int)>();
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

        auto dfs2 = function<void(int, int, ll, int)>();
        dfs2 = [&](int v, int par, ll dst, int cur_server) -> void {
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

    ~FasterKServers() {
        delete movingInstance;
    }

};


struct BinaryLifting {
    Tree tree;
    vector<vector<int>> up;
    vector<vector<ll>> dist;
    int n, LOG;

    int get_height(int v, int par = -1) {
        int ret = 0;
        for (int id : tree.g[v]) {
            int to = tree.edges[id].to;
            if (to == par)
                continue;
            ret = max(ret, 1 + get_height(to, v));
        }
        return ret;
    }

    void pre_dfs(int v, int par = -1) {
        for (int i = 1; i < LOG; ++i) {
            up[v][i] = up[up[v][i - 1]][i - 1];
            dist[v][i] = dist[v][i - 1] + dist[up[v][i - 1]][i - 1];
        }
        for (int id : tree.g[v]) {
            int to = tree.edges[id].to;
            if (to == par)
                continue;
            up[to][0] = v;
            dist[to][0] = tree.edges[id].weight;
            pre_dfs(to, v);
        }
    }

    int jump(int v, ll dst) {
        for (int i = LOG - 1; i >= 0; --i) {
            if (dist[v][i] <= dst) {
                dst -= dist[v][i];
                v = up[v][i];
            }
        }
        return v;
    }

    BinaryLifting(Tree t) {
        tree = t;
        n = tree.n;
        LOG = log2(get_height(0) + 1) + 2;
        up.resize(n, vector<int>(LOG, 0));
        dist.resize(n, vector<ll>(LOG, 0));

        pre_dfs(0);
    }

    BinaryLifting() { }
};

struct FastMoving : MovingInterface {
    Tree tree, hld_tree;
    BinaryLifting lifting;
    HLD hld;
    int n;

    void build_hld_tree() {
        vector<Edge> edgs;
        for (int v = 1; v < n; ++v) {
            if (hld.high[hld.way[v]] == v) {
                int u = hld.high[hld.way[tree.parent[v]]];
                edgs.push_back({u, v, hld.get_dist(u, v)});
            }
        }
        while (edgs.size() < n)
            edgs.push_back({n, n});
        hld_tree = Tree(edgs);
        lifting = BinaryLifting(hld_tree);
    }

    FastMoving(Tree t) {
        tree = Tree(t);
        hld = HLD(tree);
        n = tree.n;
        build_hld_tree();
    }

    int move(int from, int to, ll dst) {
        int par = hld.lca(from, to);
        int d = hld.get_dist(from, par);
        if (d < dst) {
            dst -= d;
            dst = hld.get_dist(to, par) - dst;
            from = to;
        }
        int v = from;
        d = hld.get_dist(hld.high[hld.way[v]], v);
        if (d <= dst) {
            dst -= d;
            v = hld.high[hld.way[v]];
            int u = lifting.jump(v, dst);
            dst -= hld.get_dist(v, u);
            v = u;
            if (dst != 0 && v != 0) {
                v = tree.parent[v];
                dst--;
            }
        }

        v = hld.anti_ind[hld.lf[hld.way[v]] + hld.ind[v] + dst];
        return v;

    }

    ll get_dist(int u, int v) {
        return hld.get_dist(u, v);
    }

    int lca(int u, int v) {
        return hld.lca(u, v);
    }
};

struct FastestKServers : FasterKServers{
    FastestKServers(vector<Edge> edgs, vector<int> positions): FasterKServers(edgs, positions) {
        FastMoving * tmp = new FastMoving(tree);
        setMovingInstance(tmp);
    }
};

void test_segment_tree() {
    int n = 10000;
    SegmentTree tree(n + 1);

    tree.upd(4, 10, 3);
    assert(tree.get_left(2, 5).first == 4);
    assert(tree.get_right(2, 5).first == 5);
    assert(tree.get_right(5, 100).first == 10);
    tree.upd(4, 4, 124);
    assert(tree.get_right(2, 4).second == 124);
    assert(tree.get_left(2, 5).second == 124);
    tree.upd(4, 9, 0);
    assert(tree.get_right(2, 5).first == 0);
    assert(tree.get_right(4, 10).first == 10);
    assert(tree.get_left(4, 10).first == 10);

    cout << "test_segment_tree completed succesfully" << endl;
}

void test_lca() {
    Tree tree({{0, 1},
               {1, 2},
               {2, 3},
               {3, 4}});
    HLD hld(tree);

    assert(hld.lca(0, 3) == 0);
    assert(hld.lca(4, 2) == 2);
    assert(hld.lca(4, 4) == 4);

    tree = Tree({{0, 2},
                 {0, 3},
                 {2, 4},
                 {2, 5},
                 {4, 8},
                 {4, 9},
                 {8, 10},
                 {9, 11},
                 {9, 12},
                 {3, 6},
                 {3, 7},
                 {6, 13},
                 {7, 14},
                 {7, 1}});

    hld = HLD(tree);

    assert(hld.lca(4, 4) == 4);
    assert(hld.lca(1, 3) == 3);
    assert(hld.lca(1, 8) == 0);
    assert(hld.lca(8, 5) == 2);
    assert(hld.lca(11, 12) == 9);
    assert(hld.lca(8, 9) == 4);
    assert(hld.lca(2, 8) == 2);

    cout << "test_lca completed succesfully" << endl;
}

void test_coloring() {
    Tree tree({{0, 2},
               {0, 3},
               {2, 4},
               {2, 5},
               {4, 8},
               {4, 9},
               {8, 10},
               {9, 11},
               {9, 12},
               {3, 6},
               {3, 7},
               {6, 13},
               {7, 14},
               {7, 1}});
    ColoringHLD hld(tree);

    hld.color(10, 1, 2);
    assert(hld.find_color(12, 5) == make_pair(4, 2));
    hld.color(4, 4, 0);
    assert(hld.find_color(12, 5) == make_pair(2, 2));
    assert(hld.find_color(0, 0) == make_pair(0, 2));
    hld.color(10, 1, 0);
    hld.color(8, 9, 120);
    assert(hld.find_color(5, 10) == make_pair(4, 120));
    assert(hld.find_color(5, 4) == make_pair(4, 120));
    hld.color(4, 10, 69);
    assert(hld.find_color(5, 10) == make_pair(4, 69));
    assert(hld.find_color(5, 4) == make_pair(4, 69));

    hld.color(14, 14, 111);
    assert(hld.find_color(1, 14) == make_pair(14, 111));
    assert(hld.find_color(14, 1) == make_pair(14, 111));

    cout << "test_coloring completed succesfully" << endl;
}

void test_moving() {
    Tree tree({{0, 2,  1},
               {0, 3,  1},
               {2, 4,  1},
               {2, 5,  1},
               {4, 8,  1},
               {4, 9,  1},
               {8, 10, 1},
               {9, 11, 1},
               {9, 12, 1},
               {3, 6,  1},
               {3, 7,  1},
               {6, 13, 1},
               {7, 14, 1},
               {7, 1,  1}});
    FastMoving * kek = new FastMoving(tree);
    MovingInterface * hld = kek;
    assert(hld->move(6, 0, 1) == 3);

    assert(hld->move(0, 1, 1) == 3);
    assert(hld->move(5, 11, 3) == 9);
    assert(hld->move(8, 6, 2) == 2);
    assert(hld->move(14, 14, 0) == 14);

    tree = Tree({{0, 2,  2},
                 {0, 3,  2},
                 {2, 4,  2},
                 {2, 5,  2},
                 {4, 8,  2},
                 {4, 9,  2},
                 {8, 10, 2},
                 {9, 11, 2},
                 {9, 12, 2},
                 {3, 6,  2},
                 {3, 7,  2},
                 {6, 13, 2},
                 {7, 14, 2},
                 {7, 1,  2}});
    MovingHLD * kek2 = new MovingHLD(tree);
    hld = kek2;

    assert(hld->move(8, 12, 4) == 9);
    assert(hld->move(8, 12, 3) == 4);
    assert(hld->move(2, 10, 5) == 8);
    assert(hld->move(12, 14, 9) == 0);
    assert(hld->move(12, 14, 11) == 3);
    assert(hld->move(12, 14, 13) == 7);


    cout << "test_moving completed succesfully" << endl;
}

void test_kservers() {

    KServers kServers({{0, 2,  1},
                       {0, 3,  1},
                       {2, 4,  1},
                       {2, 5,  1},
                       {4, 8,  1},
                       {4, 9,  1},
                       {8, 10, 1},
                       {9, 11, 1},
                       {9, 12, 1},
                       {3, 6,  1},
                       {3, 7,  1},
                       {6, 13, 1},
                       {7, 14, 1},
                       {7, 1,  1}},
                      {0});

    for (int i = 0; i < 15; ++i) {
        kServers.serve(i);
        assert(kServers.positions[0] == i);
    }

    kServers = KServers({{0, 2,  1},
                         {0, 3,  1},
                         {2, 4,  1},
                         {2, 5,  1},
                         {4, 8,  1},
                         {4, 9,  1},
                         {8, 10, 1},
                         {9, 11, 1},
                         {9, 12, 1},
                         {3, 6,  1},
                         {3, 7,  1},
                         {6, 13, 1},
                         {7, 14, 1},
                         {7, 1,  1}},
                        {13, 14});

    assert(kServers.serve(7) == 1);
    assert(kServers.positions[1] == 7);
    assert(kServers.serve(14) == 1);
    assert(kServers.serve(3) == 0);
    watch(kServers.positions[1]);
    assert(kServers.positions[0] == 3);
    assert(kServers.positions[1] == 7);
    cout << "test_kservers completed succesfully" << endl;
}

void test_naive() {

    Naive naive({{0, 2,  1},
                 {0, 3,  1},
                 {2, 4,  1},
                 {2, 5,  1},
                 {4, 8,  1},
                 {4, 9,  1},
                 {8, 10, 1},
                 {9, 11, 1},
                 {9, 12, 1},
                 {3, 6,  1},
                 {3, 7,  1},
                 {6, 13, 1},
                 {7, 14, 1},
                 {7, 1,  1}},
                {0});

    for (int i = 0; i < 15; ++i) {
        naive.serve(i);
        assert(naive.positions[0] == i);
    }

    naive = Naive({{0, 2,  1},
                   {0, 3,  1},
                   {2, 4,  1},
                   {2, 5,  1},
                   {4, 8,  1},
                   {4, 9,  1},
                   {8, 10, 1},
                   {9, 11, 1},
                   {9, 12, 1},
                   {3, 6,  1},
                   {3, 7,  1},
                   {6, 13, 1},
                   {7, 14, 1},
                   {7, 1,  1}},
                  {13, 14});

    assert(naive.serve(7) == 1);
    assert(naive.positions[1] == 7);
    assert(naive.serve(14) == 1);
    assert(naive.serve(3) == 0);
    watch(naive.positions[1]);
    assert(naive.positions[0] == 3);
    assert(naive.positions[1] == 7);
    cout << "test_naive completed succesfully" << endl;
}

void test_faster() {

    FasterKServers * fasterKServers = new FasterKServers({{0, 2,  1},
                                   {0, 3,  1},
                                   {2, 4,  1},
                                   {2, 5,  1},
                                   {4, 8,  1},
                                   {4, 9,  1},
                                   {8, 10, 1},
                                   {9, 11, 1},
                                   {9, 12, 1},
                                   {3, 6,  1},
                                   {3, 7,  1},
                                   {6, 13, 1},
                                   {7, 14, 1},
                                   {7, 1,  1}},
                                  {0});

    for (int i = 0; i < 15; ++i) {
        fasterKServers->serve(i);
        assert(fasterKServers->positions[0] == i);
    }
    delete fasterKServers;
    fasterKServers = new FasterKServers({{0, 2,  1},
                                     {0, 3,  1},
                                     {2, 4,  1},
                                     {2, 5,  1},
                                     {4, 8,  1},
                                     {4, 9,  1},
                                     {8, 10, 1},
                                     {9, 11, 1},
                                     {9, 12, 1},
                                     {3, 6,  1},
                                     {3, 7,  1},
                                     {6, 13, 1},
                                     {7, 14, 1},
                                     {7, 1,  1}},
                                    {13, 14});

    assert(fasterKServers->serve(7) == 1);
    assert(fasterKServers->positions[1] == 7);
    assert(fasterKServers->serve(14) == 1);
    assert(fasterKServers->serve(3) == 0);
    watch(fasterKServers->positions[1]);
    assert(fasterKServers->positions[0] == 3);
    assert(fasterKServers->positions[1] == 7);
    cout << "test_faster completed succesfully" << endl;
    delete fasterKServers;
}

void sample_testing() {
    test_segment_tree();
    test_lca();
    test_coloring();
    test_moving();
    test_kservers();
    test_naive();
    test_faster();
}

void stress_testing() {
    srand(time(0));
    int cnt = -1;
    while (true) {
        if (++cnt % 1 == 0) {
            watch(cnt);
        }
        int n = 2 + rand() % 100;
        vector<Edge> edgs;
        for (int i = 1; i < n; ++i) {
            int par = rand() % i;
            edgs.push_back({par, i, 1 + rand() % 1}); // usually 100
        }
        int k = 2 + rand() % (n - 1);
        vector<int> servers(n);
        for (int i = 0; i < n; ++i) {
            servers[i] = i;
        }
        random_shuffle(servers.begin(), servers.end());
        servers.resize(k);

        FastestKServers fastestKServers(edgs, servers);
        KServers kServers(edgs, servers);

        int q = 1 + rand() % 100;
        for (int it = 0; it < q; ++it) {
//            watch(cnt);
//            watch(it);
            int query = rand() % n;
            assert(fastestKServers.serve(query) == kServers.serve(query));
            assert(fastestKServers.positions == kServers.positions);
        }
    }
}
// TODO make O(1) LCA
int main() {
//    sample_testing();
    stress_testing();
}

