
template<class Info, class Tag>
struct LazySegmentTree {
    int n;
    vector<Info> info;
    vector<Tag> tag;
    LazySegmentTree() {

    }
    void work(int n) {
        this->n = n;
        info.resize(4 * n);
        tag.resize(4 * n);
    }
    template <typename T>
    void work(vector<T> init) {
        this->n = (int)(init.size());
        info.resize(4 * n);
        tag.resize(4 * n);
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
    LazySegmentTree(int n) : n(n), info(4 * n), tag(4 * n) {}
    template <typename T>
    LazySegmentTree(vector<T> init) : LazySegmentTree(init.size()) {
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
    void apply(int id, const Tag &v) {
        info[id] += v;
        tag[id] += v;
    }
    void push(int id) {
        apply(id << 1, tag[id]);
        apply(id << 1 | 1, tag[id]);
        tag[id] = Tag();
    }
    void modify(int id, int l, int r, int p, const Info &v) {
        if (r == l) {
            info[id] = v;
            return;
        }
        int m = l + r >> 1;
        push(id);
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
        push(id);
        return rangeQuery(id << 1, l, m, ql, qr) + rangeQuery(id << 1 | 1, m + 1, r, ql, qr);
    }
    Info rangeQuery(int ql, int qr) {
        return rangeQuery(1, 0, n - 1, ql, qr);
    }
    void rangeApply(int id, int l, int r, int ql, int qr, const Tag &v) {
        if (qr < l || r < ql) {
            return;
        }
        if (ql <= l && r <= qr) {
            apply(id, v);
            return;
        }
        int m = l + r >> 1;
        push(id);
        rangeApply(id << 1, l, m, ql, qr, v);
        rangeApply(id << 1 | 1, m + 1, r, ql, qr, v);
        pull(id);
    }
    void rangeApply(int ql, int qr, const Tag &v) {
        return rangeApply(1, 0, n - 1, ql, qr, v);
    }
};


struct Tag {
};

//lazy数据的更新
Tag& operator+=(Tag &a, Tag b) {

}

struct Info {

};

Info operator+(Info a, Info b) {

}

Info& operator+=(Info &a, Tag b) {

}