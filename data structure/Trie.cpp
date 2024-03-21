
struct Node {
    // Z val = 0 ;
    array<Node*, 2> next;
};
struct Trie {
    // Node *root = new Node();
    Node *root;
    int m;
    Trie() {
        m = 30;
        root = new Node();
    }
    Trie(int m) {
        this->m = m;
        root = new Node();
    }
    void insert(int x) {
        // for(int i=)
        for (int i = m; i >= 0; i--) {
            int nx = x >> i & 1;
            if (root->next[nx] == nullptr) {
                root->next[nx] = new Node();
            }
            root = root->next[nx];
        }
    }
};