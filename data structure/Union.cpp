struct Union {
    int N;
    vector<int> p;
    vector<int> sz;

    Union(int N) {
        init(N);
    }
    void init(int N) {
        p.resize(N);
        sz.resize(N);
        iota(all(p), 0);
        fill(all(sz), 1);
    }
    int find(int x) {
        if (x != p[x])p[x] = find(p[x]);
        return p[x];
    }

    bool merger(int from, int to) {
        from = find(from);
        to = find(to);
        if (from == to) {
            return false;
        }
        p[from] = to;
        sz[to] += sz[from];
        return true;
    }

    int getsz(int x) {
        return sz[find(x)];
    }
};