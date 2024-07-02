[toc]

## 前言

1. **注意inf，mod的初始化**
2. **正难则反**
3. **切题认真考虑一下**
4. **有一点迹象会爆ll，就得记得开ll**
5. **5.TLE，MLE算一下复杂度，复杂度对的情况下，就要进行常数优化(vector 替换掉普通数组，初始化一些重复计算(pow(2)等，组合数)**
6. **当成力扣周赛打就行**
7. **能优化尽量优化，常见的比如单调队列dp,千万别用线段树，multiset加log**
8. **dp没多大把握写递推，先尝试记忆化搜索交一发，T了再说**
9. **wa的神志不清，及时换题，解决其他题再说**
10. **可能过了很多题，心理上不要摆**

### 常见的数学公式
1. 等差数列求和
    通项公式 
    $a_1*n+d*(n*(n-1))/2$
    常见等差为1的公式
     $n*(n+1)/2$
2. 等比数列求和
    通项公式
    $a_1*(1-q^n)/(1-q)$

3. $\sum_{i=1}^{n} i=n*(n+1)/2$
4. $\sum_{i=1}^{n} i^2=n*(n+1)*(2*n+1)/6$
5. $\sum_{i=1}^{n} i^3=(n*(n+1)/2)^2$
6. 范德蒙德卷积公式 
    $\sum_{i=0}^{k} C_n^iC_{m}^{k-i} = C_{n+m}^{k} $
    推论1: $\sum_{i=-r}^{s} C_n^{r+i} C_m^{s-i}=C_{n+m}^{r+s}$

    推论2: $\sum_{i=1} ^n C _n^i C_n^{i-1} =C_{2n}^{n-1} $
    推论3：$\sum_{i=0}^n (C_n^i)^2=C_{2n}^n $
    推论4：$\sum_{i=0}^mC_n^iC_m^i=C_{n+m}^m $
7. 求切比雪夫距离 $\max(|x_i-x_j|,|y_i-y_j|)$ 等价于将点$(x,y)$转化为$(x+y,x-y)$的曼哈顿距离除2

#### 组合数
常见的就球盒模型
1.将r个无区别的球放入n个有区别的盒中，盒内数目无限制，有多少种情况？
隔板法:$C_{n+r-1}^{r}$ 

2.将r个有区别的球放入n个有区别的盒中，没有一个盒子为空，有多少种情况？
$\sum_{i=0} ^{n}(-1)^{i}*C_{n}^{i}*(n-i)^r $

3.将r个无区别的球放入n个有区别的盒中，没有一个盒子为空，有多少种情况？
$C_{r-1}^{n-1}$

4.将r个有区别的球放入n个无区别的盒中，没有一个盒子为空，有多少种情况？
$\frac{1}{n!}\sum_{i=0}^{n}(-1)^iC_n^{i}(n-i)^r $


5.将r个有区别的球放入n个有区别的盒中，盒内数目不限制，有多少种情况？
$n^r$


## 常用的一些技巧
1. 子集枚举
~~~C++
for (int i = 0; i < (1 << M); i++) {
        int u = i;
        while (true) {
            if (u == 0)
                break;
            u = (u - 1)&i;
        }
}
~~~

2. 当要进行操作要对数组进行排序
    可以考虑从**逆序对**的方面考虑

3. 树上路径问题 
    从LCA，树链剖分入手

4. 贪心大胆发挥想象空间，有技巧的暴力，分类讨论等等

5. 图上问题可以考虑转化成树考虑
6. 子数组问题常常固定左端点，枚举右端点
   

## 数据结构

1. lazy线段树根据参考，赛时可简写

~~~C++

    template<class Info, class Tag>
    struct LazySegmentTree {
        int n;
        vector<Info> info;
        vector<Tag> tag;
        LazySegmentTree() {

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
~~~
2. ST表常见用于静态查询($gcd,lcm,或，与，LCA ,max,min$)
~~~C++
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
};
~~~
3.树状数组 (1. 区间查询，单点修改。2，区间修，单点查)
 ~~~C++ 
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
~~~

4. 树链剖分 看着是树剖，但也可以不用写,除非有十足的把握
   
~~~C++
// in 表示重新编号的ID
// seq[x] 以x为id的节点是什么
// top[x] 表示这个链的链头
struct HLD {
    int n;
    std::vector<int> siz, top, dep, parent, in, out, seq;
    std::vector<std::vector<int>> adj;
    int cur;

    HLD() {}

    HLD(int n) {
        init(n);
    }

    void init(int n) {
        this->n = n;
        siz.resize(n);
        top.resize(n);
        dep.resize(n);
        parent.resize(n);
        in.resize(n);
        out.resize(n);
        seq.resize(n);
        cur = 0;
        adj.assign(n, {});
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    void work(int root = 0) {
        top[root] = root;
        dep[root] = 0;
        parent[root] = -1;
        dfs1(root);
        dfs2(root);
    }

    void dfs1(int u) {
        if (parent[u] != -1) {
            adj[u].erase(std::find(adj[u].begin(), adj[u].end(), parent[u]));
        }

        siz[u] = 1;
        for (auto &v : adj[u]) {
            parent[v] = u;
            dep[v] = dep[u] + 1;
            dfs1(v);
            siz[u] += siz[v];
            if (siz[v] > siz[adj[u][0]]) {
                std::swap(v, adj[u][0]);
            }
        }
    }

    void dfs2(int u) {
        in[u] = cur++;
        seq[in[u]] = u;
        for (auto v : adj[u]) {
            top[v] = v == adj[u][0] ? top[u] : v;
            dfs2(v);
        }
        out[u] = cur;
    }

    // 求LCA
    int lca(int u, int v) {
        while (top[u] != top[v]) {
            if (dep[top[u]] > dep[top[v]]) {
                u = parent[top[u]];
            } else {
                v = parent[top[v]];
            }
        }
        return dep[u] < dep[v] ? u : v;
    }

    int dist(int u, int v) {
        return dep[u] + dep[v] - 2 * dep[lca(u, v)];
    }

    // 向上跳K步
    int jump(int u, int k) {
        if (dep[u] < k) {
            return -1;
        }

        int d = dep[u] - k;

        while (dep[top[u]] > d) {
            u = parent[top[u]];
        }

        return seq[in[u] - dep[u] + d];
    }

    bool isAncester(int u, int v) {
        return in[u] <= in[v] && in[v] < out[u];
    }
    //返回 u->v 路径上 v的前一个节点
    int rootedParent(int u, int v) {
        std::swap(u, v);
        if (u == v) {
            return u;
        }
        if (!isAncester(u, v)) {
            return parent[u];
        }
        auto it = std::upper_bound(adj[u].begin(), adj[u].end(), v, [&](int x, int y) {
            return in[x] < in[y];
        }) - 1;
        return *it;
    }

    int rootedSize(int u, int v) {
        if (u == v) {
            return n;
        }
        if (!isAncester(v, u)) {
            return siz[v];
        }
        return n - siz[rootedParent(u, v)];
    }

    int rootedLca(int a, int b, int c) {
        return lca(a, b) ^ lca(b, c) ^ lca(c, a);
    }
};
~~~

5. 主席树求K小
~~~C++
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
~~~

## 字符串相关模板
1. 字符串哈希
**注意取模**
~~~C++
struct Hash_String
{
    // string s;
    vector<Z> h, p;
    Hash_String() {

    }
    Hash_String(string t) {
        work(t);
    }
    void work(string t) {
        t = " " + t;
        int N = len(t);
        h.resize(N + 1);
        p.resize(N + 1);
        p[0] = 1;
        for (int i = 1; i < N; i ++ )
        {
            h[i] = h[i - 1] * B + t[i];
            p[i] = p[i - 1] * B;
        }
    }
    Z get(int l, int r) {
        return h[r] - h[l - 1] * p[r - l + 1];
    }
};
~~~
2. KMP算法
~~~C++
void computeLPSArray( std::string& pattern, std::vector<int>& lps) {
    int len = 0;
    lps[0] = 0;
    int i = 1;
    while (i < pattern.length()) {
        if (pattern[i] == pattern[len]) {
            len++;
            lps[i] = len;
            i++;
        } else {
            if (len != 0) {
                len = lps[len - 1];
            } else {
                lps[i] = 0;
                i++;
            }
        }
    }
}

// 字符串匹配求模式项与匹配项
vector<int>  KMPSearch( std::string& pattern,  std::string& text) {
    int M = pattern.length();
    int N = text.length();

    std::vector<int> lps(M);
    computeLPSArray(pattern, lps);
    vector<int> ans;
    int i = 0;
    int j = 0;
    while (i < N) {
        if (pattern[j] == text[i]) {
            j++;
            i++;
        }
        if (j == M) {
            ans.push_back(i-j);
            j = lps[j - 1];
        } else if (i < N && pattern[j] != text[i]) {
            if (j != 0) {
                j = lps[j - 1];
            } else {
                i = i + 1;
            }
        }
    }
    return ans;
}
~~~

3. Z函数求后缀与整个字符串的LCP
~~~C++
vector<int> z_function_trivial(string s) {
    int n = (int)s.length();
    vector<int> z(n);
    for (int i = 1; i < n; ++i)
        while (i + z[i] < n && s[z[i]] == s[i + z[i]])
            ++z[i];
    return z;
}
~~~

## 数论
没啥好说的，真的没啥好说的，听天由命吧
1. exgcd
~~~C++
int gcd(int a, int b, int& x, int& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    int x1, y1;
    int d = gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - a / b * y1;
    return d;
}

bool find_solution(int a, int b, int c, int& x0, int& y0) {
    int x, y;
    int g = gcd(abs(a), abs(b), x, y);

    if (c % g != 0) {
        return false;
    }

    x0 = x * c / g;
    y0 = y * c / g;

    if (a < 0) x0 = -x0;
    if (b < 0) y0 = -y0;

    return true;
}
~~~

## 图论
1. 强调一下最短路算法，很多情况都可以考虑最短路算法(dijkstra) ,并且明确是个最短路算法的时候(dijkstra)，一定要发挥想象空间，根据题意进行魔改，分层图

2. 缩点(算强连通分量)
~~~C++
struct Scc
{
    // g是新图
    // scc  原先的点现在所处的位置
    // sz 表示现在的点的联通大小
    vector<vector<int>> g;
    vector<int> dfn, low, stk, scc, sz;
    vector<bool> instak;
    
    int tot = 0, top = 0, cnt = 0;
    Scc(int N, vector<vector<int>> &e) {
        dfn.resize(N + 1);
        scc.resize(N + 1);
        sz.resize(N + 1);
        instak.resize(N + 1);
        stk.resize(N + 1);
        low.resize(N + 1);

        auto dfs = [&](auto dfs, int x)->void{
            dfn[x] = low[x] = ++tot;
            stk[++top] = x;
            instak[x] = true;
            for (int y : e[x]) {
                if (dfn[y] == 0) {
                    dfs(dfs, y);
                    low[x] = min(low[x], low[y]);
                } else if (instak[y]) {
                    low[x] = min(low[x], dfn[y]);
                }
            }
            if (dfn[x] == low[x]) {
                int y;
                ++cnt;
                do {
                    y = stk[top--];
                    instak[y] = false;
                    scc[y] = cnt;
                    ++sz[cnt];
                } while (y != x);
            }
        };
        for (int i = 1; i <= N; i++) {
            if (dfn[i] == 0) {
                dfs(dfs, i);
            }
        }

        g.resize(cnt + 1);
        for (int i = 1; i <= N; i++) {
            for (int x : e[i]) {
                if (scc[i] != scc[x]) {
                    g[scc[i]].push_back(scc[x]);
                }
            }
        }
    }
};
~~~
1. 并查集
~~~C++
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
~~~