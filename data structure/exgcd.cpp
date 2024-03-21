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