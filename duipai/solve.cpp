#include<bits/stdc++.h>

using namespace std;

#define endl '\n'
#define  len(x) ((int)((x).size()))
#define PQ priority_queue
#define all(x) (x).begin(),(x).end()

using ll = long long;
using ld = long double;

mt19937_64 mrand(random_device{}());

const array<int, 8> dx{0, -1, 0, 1, 1, 1, -1, -1};
const array<int, 8> dy{1, 0, -1, 0, 1, -1, 1, -1};


signed main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);
    int a, b;
    cin >> a >> b;
    
    cout << a + b << endl;


    return 0;
}
