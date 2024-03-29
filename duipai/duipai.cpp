#include <bits/stdc++.h>

using namespace std;

int main() {
    int t = 0;
    while (true) {
        cout << "test£º" << t++ << endl;
        system("data.exe > data.in");
        system("std.exe < data.in > std.out");
        system("java -jar myCode.jar < data.in > solve.out");
        if (system("fc std.out solve.out > diff.log")) {
            cout << "WA\n";
            break;
        }
        cout << "AC\n";
    }
}
