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
            // std::cout << "Pattern found at index " << i - j << std::endl;
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