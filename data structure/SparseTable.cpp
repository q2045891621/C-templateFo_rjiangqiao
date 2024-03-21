template <typename T>
struct SparseTable {
    vector<vector<T>> st;
    vector<int> logT;
    using func = function<T(T, T)>;
    func f;
    int N;
    // 传递数组，函数
    SparseTable(vector<T> &a, func f1) {
        N = len(a);
        int logN = log2(N) + 1;
        f = f1;
        st.resize(N);
        logT.resize(N + 1);
        for (int i = 0; i < N; i++) {
            st[i].resize(logN);
            st[i][0] = a[i];
        }
        for (int i = 2; i <= N; i++) {
            logT[i] = logT[i >> 1] + 1;
        }
        for (int j = 1; (1 << j) <= N; j++) {
            for (int i = 0; i + (1 << j) - 1 < N; i++) {
                st[i][j] = f(st[i][j - 1], st[i + (1 << (j - 1))][j - 1]);
            }
        }
    }

    // 区间查询
    T query(int left, int right) {
        // left<=right
        assert(left <= right);
        int k = logT[right - left + 1];
        return f(st[left][k], st[right - (1 << k) + 1][k]);
    }

    int get_left(int u) {
        int left = 0;
        int right = u;
        T cur = query(u, u);
        while (left < right) {
            int mid = (left + right) / 2;
            if (query(mid, u) == cur) {
                right = mid;
            } else left = mid + 1;
        }
        return left;
    }

    int get_right(int u) {
        int left = u;
        int right = N - 1 ;
        T cur = query(u, u);
        while (left < right) {
            int mid = (1 + left + right) / 2;
            if (query(u, mid) == cur) {
                left = mid;
            } else right = mid - 1;
        }
        return left;
    }
};