// 区间K小,下标从1开始
struct segment_tree_node
{
    int lc, rc, sum;
};
struct Presidents_tree {
    vector<segment_tree_node> t;
    vector<int> next;
    vector<int> rt;
    int cnt = 0;
    int sz;

    Presidents_tree(vector<int> &a) {
        next = a;
        int N = (int)(a.size());
        sort(all(next));
        t.resize(N * 30);
        rt.resize(N + 100);
        sort(next.begin() + 1, next.end());
        sz = unique(next.begin() + 1, next.end()) - next.begin() - 1;
        for (int i = 1; i < N; i++) {
            int x = lower_bound(next.begin() + 1, next.begin() + sz, a[i]) - next.begin();
            rt[i] = change(rt[i - 1], 1, sz, x);
        }
    }
    int change(int k, int l, int r, int x) {
        if (x > r || x < l)
            return 0;
        int p = ++cnt, mid = (l + r) >> 1;
        t[p] = t[k];
        if (l == r) {
            t[p].sum++;
            return p;
        }
        if (x <= mid)
            t[p].lc = change(t[k].lc, l, mid, x);
        else
            t[p].rc = change(t[k].rc, mid + 1, r, x);
        t[p].sum = t[t[p].lc].sum + t[t[p].rc].sum;
        return p;
    }
    int query(int k1, int k2, int l, int r, int x) {
        int mid = (l + r) >> 1;
        if (l == r)
            return next[l];
        int ans = 0;
        if (t[t[k2].lc].sum - t[t[k1].lc].sum >= x)
            return query(t[k1].lc, t[k2].lc, l, mid, x);
        else
            return query(t[k1].rc, t[k2].rc, mid + 1, r, x - (t[t[k2].lc].sum - t[t[k1].lc].sum));
    }

    int query(int L, int R, int K) {
        return query(rt[L - 1], rt[R], 1, sz, K);
    }

};