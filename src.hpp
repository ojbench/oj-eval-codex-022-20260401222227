#include <bits/stdc++.h>
using namespace std;

// The OJ provides this function; we only declare it.
int query(int x, int y, int z);

static const int MOD = 998244353;
static const int MUL = 233;

// Helper to compute hash of recovered A
static int compute_hash(const vector<long long>& A) {
    long long ret = 0;
    for (int i = (int)A.size() - 1; i >= 0; --i) {
        ret = (ret + (A[i] % MOD + MOD) % MOD) % MOD;
        ret = (ret * MUL) % MOD;
    }
    return (int)ret;
}

// Try to find indices of global min and max
static pair<int,int> find_extremes(int n) {
    if (n <= 2) return {1, n};
    int a = 1, b = 2;
    // First pass: use pair (1,2)
    long long mnv = LLONG_MAX, mxv = LLONG_MIN; int mnidx = -1, mxidx = -1;
    for (int i = 3; i <= n; ++i) {
        int s = query(a, b, i);
        if ((long long)s < mnv) { mnv = s; mnidx = i; }
        if ((long long)s > mxv) { mxv = s; mxidx = i; }
    }
    a = mnidx; b = mxidx;
    if (a == -1 || b == -1) {
        // fallback: pick any different indices
        a = 1; b = 2;
    }
    // Second pass: refine
    mnv = LLONG_MAX; mxv = LLONG_MIN; mnidx = -1; mxidx = -1;
    for (int i = 1; i <= n; ++i) if (i != a && i != b) {
        int s = query(a, b, i);
        if ((long long)s < mnv) { mnv = s; mnidx = i; }
        if ((long long)s > mxv) { mxv = s; mxidx = i; }
    }
    int a2 = mnidx, b2 = mxidx;
    if (a2 == -1 || b2 == -1) return {a, b};
    // Third pass: final candidates
    a = a2; b = b2;
    return {a, b};
}

int guess(int n, int Taskid) {
    vector<long long> A(n + 1, 0);
    auto ext = find_extremes(n);
    int imin = ext.first, imax = ext.second;

    // Choose an index between extremes to get S = A_min + A_max
    int t_between = -1;
    for (int i = 1; i <= n; ++i) if (i != imin && i != imax) { t_between = i; break; }
    if (t_between == -1) t_between = (imin == 1 ? 2 : 1);
    int S = query(imin, imax, t_between);

    // Select three indices p, q, r (attempt non-extremes when possible)
    vector<int> nonext;
    for (int i = 1; i <= n; ++i) if (i != imin && i != imax) nonext.push_back(i);
    int p = -1, q = -1, r = -1;
    if ((int)nonext.size() >= 3) { p = nonext[0]; q = nonext[1]; r = nonext[2]; }
    else if ((int)nonext.size() == 2) { p = nonext[0]; q = nonext[1]; r = imin; }
    else { p = imin; q = imax; r = t_between; }

    // Using extremes, compute pairwise sums among p, q, r
    auto sum2 = [&](int x, int y) -> long long {
        return (long long)query(imin, x, y) + (long long)query(imax, x, y) - (long long)S; // Ax + Ay
    };

    long long spq = sum2(p, q);
    long long spr = sum2(p, r);
    long long sqr = sum2(q, r);

    long long Ap = (spq + spr - sqr) / 2;
    long long Aq = (spq + sqr - spr) / 2;
    long long Ar = (spr + sqr - spq) / 2;

    A[p] = Ap; A[q] = Aq; A[r] = Ar;

    // Recover others: for each j, (imin,j,p)+(imax,j,p)-S = Aj + Ap
    for (int j = 1; j <= n; ++j) {
        if (j == imin || j == imax || j == p || j == q || j == r) continue;
        long long u = (long long)query(imin, j, p) + (long long)query(imax, j, p) - (long long)S;
        A[j] = u - Ap;
    }

    // Recover A_min and A_max using p, q
    long long mxpq = max(Ap, Aq), mnpq = min(Ap, Aq);
    A[imin] = (long long)query(imin, p, q) - mxpq;
    A[imax] = (long long)query(imax, p, q) - mnpq;

    // Compute hash
    return compute_hash(A);
}
