struct BipMatch
{
    // vector<vector<int>> g;
    //N表示男生数量，M表示女生数量

    // int N, M;
    int ans = 0;
    BipMatch(vector<vector<int>> &g, int M) {
        vector <bool> vis(M);
        int N = len(g);
        vector<int> match(M, -1);
        auto dfs = [&](auto dfs, int u)->bool{
            for (int ne : g[u]) {
                if (vis[ne])continue;
                vis[ne] = 1;
                if (match[ne] == -1 || dfs(dfs, match[ne])) {
                    match[ne] = u;
                    return 1;
                }
            }
            return 0;
        };

        for (int i = 0; i < N; i++) {
            fill(all(vis), 0);
            if (dfs(dfs, i)) {
                ans++;
            }
        }
    }
    int get() {
        return ans;
    }
};