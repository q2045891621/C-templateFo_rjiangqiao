template<class Info>
struct SegmentTree {
    int n;
    vector<Info> info;
    SegmentTree(int n) : n(n), info(4 * n) {}
    template <typename T>
    SegmentTree(vector<T> init) : SegmentTree(init.size()) {
        auto build = [&](auto build, int id, int l, int r) -> void {
            if (r == l) {
                info[id] = init[l];
                return;
            }
            int m = l + r >> 1;
            build(build, id << 1, l, m);
            build(build, id << 1 | 1, m + 1, r);
            pull(id);
        };
        build(build, 1, 0, n - 1);
    }
    void pull(int id) {
        info[id] = info[id << 1] + info[id << 1 | 1];
    }
    void modify(int id, int l, int r, int p, const Info &v) {
        if (r == l) {
            info[id] += v;
            return;
        }
        int m = l + r >> 1;
        if (p <= m) {
            modify(id << 1, l, m, p, v);
        } else {
            modify(id << 1 | 1, m + 1, r, p, v);
        }
        pull(id);
    }
    void modify(int p, const Info &v) {
        modify(1, 0, n - 1, p, v);
    }
    Info rangeQuery(int id, int l, int r, int ql, int qr) {
        if (qr < l || r < ql) {
            return Info();
        }
        if (ql <= l && r <= qr) {
            return info[id];
        }
        int m = l + r >> 1;

        return rangeQuery(id << 1, l, m, ql, qr) + rangeQuery(id << 1 | 1, m + 1, r, ql, qr);
    }
    Info rangeQuery(int ql, int qr) {
        return rangeQuery(1, 0, n - 1, ql, qr);
    }

};

//存放节点数据,注意初始化
struct Info {

};
// 相当于两个区间的合并
Info operator+(Info a, Info b) {

}
// 单点修改的逻辑
Info &operator+=(Info &a, Info b) {

}