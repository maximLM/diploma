#include <bits/stdc++.h>

#include "includes/Naive.h"
#include "includes/FirstSolution.h"
#include "includes/SecondSolution.h"
#include "includes/ThirdSolution.h"
#include "includes/structures/FastMoving.h"

#define watch(x) cout << (#x) << " = " << (x) << endl;
// comment
using namespace std;


typedef long long ll;
typedef pair<int, int> pii;

const ll oo = 1e9 + 10;


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

    FirstSolution kServers({{0, 2,  1},
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

    kServers = FirstSolution({{0, 2,  1},
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

    SecondSolution * fasterKServers = new SecondSolution({{0, 2,  1},
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
    fasterKServers = new SecondSolution({{0, 2,  1},
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
    int cnt = -1;
    while (true) {
        if (++cnt % 100 == 0) {
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

        ThirdSolution fastestKServers(edgs, servers);
        Naive naive(edgs, servers);

        int q = 1 + rand() % 100;
        for (int it = 0; it < q; ++it) {
//            watch(cnt);
//            watch(it);
            int query = rand() % n;
            assert(fastestKServers.serve(query) == naive.serve(query));
            assert(fastestKServers.positions == naive.positions);
        }
    }
}

void benchmark() {
    srand(time(0));
    int n = 1000;
    freopen("output.txt", "w", stdout);

    for (int k = 20; k <= 1000; k += 10) {
        cerr << int((double)k / 1000 * 100) << "%" << endl;
        vector<Edge> edgs;
        for (int i = 1; i < n; ++i) {
            int par = rand() % i;
            edgs.push_back({par, i, 1 + rand() % 1}); // usually 100
        }
        vector<int> servers(n);
        for (int i = 0; i < n; ++i) {
            servers[i] = i;
        }
        random_shuffle(servers.begin(), servers.end());
        servers.resize(k);

        Naive naive(edgs, servers);
        FirstSolution kservers(edgs, servers);
        SecondSolution fasterKServers(edgs, servers);
        ThirdSolution fastestKServers(edgs, servers);

        int q = 10;
        vector<double> finalTimes(4);
        for (int it = 0; it < q; ++it) {
            int v = rand() % n;
            vector<double> tms;
            clock_t start;
            start = clock();
            naive.serve(v);
            tms.push_back(double(clock() - start) / CLOCKS_PER_SEC);


            start = clock();
            kservers.serve(v);
            tms.push_back(double(clock() - start) / CLOCKS_PER_SEC);


            start = clock();
            fasterKServers.serve(v);
            tms.push_back(double(clock() - start) / CLOCKS_PER_SEC);


            start = clock();
            fastestKServers.serve(v);
            tms.push_back(double(clock() - start) / CLOCKS_PER_SEC);

            for (int i = 0; i < tms.size(); ++i) {
                finalTimes[i] += tms[i];
            }
        }
        cout << n << ' ' << k << ' ';
        for (auto & x : finalTimes) {
            x /= q;
            cout << x << ' ';
        }
        cout << '\n';
    }
}

int main() {
//    sample_testing();
    stress_testing();
//    benchmark();
}

