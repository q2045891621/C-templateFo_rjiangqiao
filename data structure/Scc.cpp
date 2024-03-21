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
