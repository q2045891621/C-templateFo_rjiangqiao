template<typename T>
struct Comb {
    int N;
    vector<vector<T>> com;
    Comb(int n) {
        N = n;
        com.resize(N + 1);
        for (int i = 0; i <= N; i++) {
            com[i].resize(N + 1);
        }
        for (int i = 0; i <= N; i++) {

            // cout<<i<<endl;
            com[i][0] = com[i][i] = 1;
            for (int j = 1; j < i; j++) {
                // cout<<i<<" "<<j<<endl;
                com[i][j] = com[i - 1][j] + com[i - 1][j - 1];
            }
        }
    }
    T get(int n, int m) {
        return com[n][m];
    }

};