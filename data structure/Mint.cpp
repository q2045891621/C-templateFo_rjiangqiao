template<typename T, T MOD = 1000000007>
struct Mint {
	static constexpr T mod = MOD;
	T v;
	Mint(): v(0) {}
	Mint(signed v): v(v) {}
	Mint(long long t) {v = t % MOD; if (v < 0) v += MOD;}

	Mint pow(long long k) {
		Mint res(1), tmp(v);
		while (k) {
			if (k & 1) res *= tmp;
			tmp *= tmp;
			k >>= 1;
		}
		return res;
	}

	static Mint add_identity() {return Mint(0);}
	static Mint mul_identity() {return Mint(1);}

	Mint inv() {return pow(MOD - 2);}

	Mint& operator+=(Mint a) {v += a.v; if (v >= MOD)v -= MOD; return *this;}
	Mint& operator-=(Mint a) {v += MOD - a.v; if (v >= MOD)v -= MOD; return *this;}
	Mint& operator*=(Mint a) {v = 1LL * v * a.v % MOD; return *this;}
	Mint& operator/=(Mint a) {return (*this) *= a.inv();}

	Mint operator+(Mint a) const {return Mint(v) += a;};
	Mint operator-(Mint a) const {return Mint(v) -= a;};
	Mint operator*(Mint a) const {return Mint(v) *= a;};
	Mint operator/(Mint a) const {return Mint(v) /= a;};

	Mint operator-() const {return v ? Mint(MOD - v) : Mint(v);}

	bool operator==(const Mint a)const {return v == a.v;}
	bool operator!=(const Mint a)const {return v != a.v;}
	bool operator <(const Mint a)const {return v < a.v;}

	// find x s.t. a^x = b
	static T log(T a, T b) {
		const T sq = 40000;
		unordered_map<T, T> dp;
		dp.reserve(sq);
		Mint res(1);
		for (ll r = 0; r < sq; r++) {
			if (!dp.count(res.v)) dp[res.v] = r;
			res *= a;
		}
		Mint p = Mint(a).inv().pow(sq);
		res = b;
		for (ll q = 0; q <= MOD / sq + 1; q++) {
			if (dp.count(res.v)) {
				T idx = q * sq + dp[res.v];
				if (idx > 0) return idx;
			}
			res *= p;
		}
		assert(0);
		return T(-1);
	}

	static Mint comb(long long n, ll k) {
		Mint num(1), dom(1);
		for (ll i = 0; i < k; i++) {
			num *= Mint(n - i);
			dom *= Mint(i + 1);
		}
		return num / dom;
	}
};
template<typename T,T MOD> constexpr T Mint<T, MOD>::mod;
template<typename T,T MOD>
ostream& operator<<(ostream &os,Mint<T, MOD> m){os<<m.v;return os;}
using Z = Mint<int>;