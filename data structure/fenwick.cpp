
template <typename T>
class fenwick {
public:
    vector<T> fenw;
    int n;
    vector<pair<int, T>> data;
    fenwick(int _n) : n(_n) {
        fenw.resize(n);
    }
    void modify(int x, T v) {
        change(x, v, 1);
    }

    T get(int x) {
        T v{};
        while (x >= 0) {
            v += fenw[x];
            x = (x & (x + 1)) - 1;
        }
        return v;
    }
    T query(int L, int R) {
        // T v{};
        return get(R) - get(L - 1);
    }
    void change(int x, T v, bool ok) {
        if (ok) {
            data.push_back({x, v});
        }
        while (x < n) {
            fenw[x] += v;
            x |= (x + 1);
        }
    }
    void clear() {
        for (auto [x, y] : data) {
            change(x, y * -1, 0);
        }
        data.clear();
    }
};