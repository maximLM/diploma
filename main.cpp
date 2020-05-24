#include <bits/stdc++.h>

#define watch(x) cout << (#x) << " = " << (x) << endl;

using namespace std;

typedef long long ll;
typedef pair<int, int> pii;

const ll oo = 1e9 + 10;

struct Edge {
    int from, to, weight;
};

struct Tree {
    int n;
    vector<Edge> edges;
    vector<vector<int>> g;
    vector<int> parent, parent_edge;

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

    Tree(vector<Edge> eds) {
        n = eds.size() + 1;
        g.resize(n);
        parent.resize(n);
        parent_edge.resize(n);
        for (auto e : eds) {
            g[e.from].push_back(edges.size());
            edges.push_back(e);
            swap(e.from, e.to);
            g[e.from].push_back(edges.size());
            edges.push_back(e);
        }
        dfs(0);
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
    SegmentTree segmentTree;
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
        segmentTree = SegmentTree(n);
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
                dst -=  tree.edges[tree.parent_edge[from]].weight;
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

    int move_down(int from, int to, ll & dst) {
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
        int l = way[v] == way[to] ? lf[way[to]] + ind[to] - 1: lf[way[v]] + ind[tree.parent[highs.back()]] - 1;
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

    int move_dist(int from, int to, ll dst) {
        int par = lca(from, to);
        int ans = move_up(from, par, dst);
        if (ans == par)
            move_down(par, to, dst);
        return ans;
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
            if (cur.first)
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
            if (cur.first)
                ret = cur;
            to = tree.parent[high[wind]];
        }

        int wind = way[to];
        pii cur = segmentTree.get_right(lf[wind] + ind[to], lf[wind] + way[from]);
        if (cur.first)
            ret = cur;
        return ret;
    }

    pii find_color(int from, int to) {
        int par = lca(from, to);
        pii ret = find_up(from, par);
        if (!ret.first)
            ret = find_down(par, to);
        assert(ret.first);
        ret.first = anti_ind[ret.first];
        return ret;
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
    HLD hld(tree);

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
    Tree tree({{0, 2, 1},
               {0, 3, 1},
               {2, 4, 1},
               {2, 5, 1},
               {4, 8, 1},
               {4, 9, 1},
               {8, 10, 1},
               {9, 11, 1},
               {9, 12, 1},
               {3, 6, 1},
               {3, 7, 1},
               {6, 13, 1},
               {7, 14, 1},
               {7, 1, 1}});
    HLD hld(tree);

}


int main() {
    test_segment_tree();
    test_lca();
    test_coloring();
}