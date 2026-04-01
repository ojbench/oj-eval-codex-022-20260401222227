#include <bits/stdc++.h>
using namespace std;

// The OJ provides this function; we only declare it.
int query(int x, int y, int z);

static const int MOD = 998244353;
static const int MUL = 233;

// Helper to compute hash of recovered A
static int compute_hash(const vector<long long>& A, int n) {
    long long ret = 0;
    for (int i = n; i >= 1; --i) {
        ret = (ret + (A[i] % MOD + MOD) % MOD) % MOD;
        ret = (ret * MUL) % MOD;
    }
    return (int)ret;
}

// Try to find indices of global min and max
static void find_extremes_and_S(int n, int &imin, int &imax, int &mid, int &S) {
    // Use 4 indices to seed extremes
    int idx[4] = {1, 2, 3, 4};
    int s_excl[4];
    s_excl[0] = query(2, 3, 4); // exclude 1
    s_excl[1] = query(1, 3, 4); // exclude 2
    s_excl[2] = query(1, 2, 4); // exclude 3
    s_excl[3] = query(1, 2, 3); // exclude 4
    int argmin = 0, argmax = 0;
    for (int i = 1; i < 4; ++i) {
        if (s_excl[i] < s_excl[argmin]) argmin = i;
        if (s_excl[i] > s_excl[argmax]) argmax = i;
    }
    // Excluding max yields minimum; excluding min yields maximum
    imax = idx[argmin];
    imin = idx[argmax];
    // pick a middle index among the remaining
    mid = -1;
    for (int i = 0; i < 4; ++i) if (idx[i] != imin && idx[i] != imax) { mid = idx[i]; break; }
    if (mid == -1) mid = (imin == 1 ? 2 : 1);
    // Compute S = Amin + Amax
    S = query(imin, imax, mid);
    // Scan through all other indices to update extremes
    for (int i = 1; i <= n; ++i) {
        if (i == imin || i == imax || i == mid) continue;
        int t = query(imin, imax, i);
        if (t < S) { imin = i; S = t; }
        else if (t > S) { imax = i; S = t; }
    }
}

int guess(int n, int Taskid) {
    vector<long long> A(n + 1, 0);
    int imin = 1, imax = 2, mid = 3, S = 0;
    find_extremes_and_S(n, imin, imax, mid, S);

    // Select three indices p, q, r (attempt non-extremes when possible)
    vector<int> nonext;
    for (int i = 1; i <= n; ++i) if (i != imin && i != imax) nonext.push_back(i);
    int p = -1, q = -1, r = -1;
    if ((int)nonext.size() >= 3) { p = nonext[0]; q = nonext[1]; r = nonext[2]; }
    else if ((int)nonext.size() == 2) { p = nonext[0]; q = nonext[1]; r = imin; }
    else { p = imin; q = imax; r = mid; }

    // Using extremes, compute pairwise sums among p, q, r
    auto sum2 = [&](int x, int y) -> long long {
        return (long long)query(imin, x, y) + (long long)query(imax, x, y) - (long long)S; // Ax + Ay
    };

    long long v_im_pq = query(imin, p, q);
    long long v_iM_pq = query(imax, p, q);
    long long spq = v_im_pq + v_iM_pq - S;
    long long v_im_pr = query(imin, p, r);
    long long v_iM_pr = query(imax, p, r);
    long long spr = v_im_pr + v_iM_pr - S;
    long long v_im_qr = query(imin, q, r);
    long long v_iM_qr = query(imax, q, r);
    long long sqr = v_im_qr + v_iM_qr - S;

    long long Ap = (spq + spr - sqr) / 2;
    long long Aq = (spq + sqr - spr) / 2;
    long long Ar = (spr + sqr - spq) / 2;

    A[p] = Ap; A[q] = Aq; A[r] = Ar;

    // Recover others using dynamic second-maximum strategy to cut queries
    // We already have Ap, Aq, Ar; choose current second-maximum among them
    int smax_idx = p; long long Asmax = Ap;
    if (Aq > Asmax) { Asmax = Aq; smax_idx = q; }
    if (Ar > Asmax) { Asmax = Ar; smax_idx = r; }

    // Compute A_min and A_max now (needed for single-query recovery)
    long long mxpq = max(Ap, Aq), mnpq = min(Ap, Aq);
    A[imin] = v_im_pq - mxpq;
    A[imax] = v_iM_pq - mnpq;

    long long Amin = A[imin], Amax = A[imax];

    for (int j = 1; j <= n; ++j) {
        if (j == imin || j == imax || j == p || j == q || j == r) continue;
        long long u = query(imax, smax_idx, j); // Amax + min(Asmax, Aj)
        if (u < Amax + Asmax) {
            A[j] = u - Amax; // Aj < Asmax
        } else {
            // Aj >= Asmax; determine Aj and update second-maximum
            long long v = query(imin, smax_idx, j); // Amin + max(Asmax, Aj) = Amin + Aj
            A[j] = v - Amin;
            if (A[j] > Asmax) { Asmax = A[j]; smax_idx = j; }
        }
    }

    // Compute hash
    return compute_hash(A, n);
}
