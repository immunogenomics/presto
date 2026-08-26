#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <thread>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;


// [[Rcpp::export]]
arma::mat cpp_sumGroups_dgc(const arma::vec& x, const arma::uvec& p, 
                            const arma::vec& i, unsigned ncol, 
                            const arma::uvec& groups, unsigned ngroups) {
    // Here, columns are genes
    arma::mat res = arma::zeros<arma::mat>(ngroups, ncol);
    for (unsigned c = 0; c < ncol; c++) {
        for (unsigned j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            // group_map gives the group num of that row
            res(groups[i[j]], c) += x[j];
        }
    }    
    return res;
}

// [[Rcpp::export]]
arma::mat cpp_sumGroups_dgc_T(const arma::vec& x, const arma::vec& p, 
                              const arma::vec& i, int ncol, int nrow, 
                              const arma::uvec& groups, int ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, nrow);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            // group_map gives the group num of that row
            res(groups[c], i[j]) += x[j];
        }
    }    
    return res;
}

// [[Rcpp::export]]
arma::mat cpp_sumGroups_dense(const arma::mat& X, const arma::uvec& groups,
                              unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
    for (unsigned r = 0; r < X.n_rows; r++) {
        res.row(groups[r]) += sum(X.row(r), 0);
    }    
    return res;
}

// [[Rcpp::export]]
arma::mat cpp_sumGroups_dense_T(const arma::mat& X, const arma::uvec& groups,
                                unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
    for (unsigned c = 0; c < X.n_cols; c++) {
        res.row(groups[c]) += sum(X.col(c), 1).t();
    }    
    return res;
}


// [[Rcpp::export]]
arma::mat cpp_nnzeroGroups_dense(const arma::mat& X, const arma::uvec& groups,
                                 unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
    for (unsigned c = 0; c < X.n_cols; c++) {
        for (unsigned r = 0; r < X.n_rows; r++) {
            if (X(r, c) != 0)
                res(groups[r], c)++;
        }
    }    
    return res;
}

// [[Rcpp::export]]
arma::mat cpp_nnzeroGroups_dense_T(const arma::mat& X,
                                   const arma::uvec& groups, 
                                   unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
    for (unsigned c = 0; c < X.n_cols; c++) {
        for (unsigned r = 0; r < X.n_rows; r++) {
            if (X(r, c) != 0)
                res(groups[c], r)++;
//             res.row(groups[c]) += sum(X.col(c), 1).t();
        }
    }    
    return res;
}



// [[Rcpp::export]]
arma::mat cpp_nnzeroGroups_dgc(const arma::uvec& p, const arma::vec& i,
                               unsigned ncol, const arma::uvec& groups,
                               unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, ncol);
    for (unsigned c = 0; c < ncol; c++) {
        for (unsigned j = p[c]; j < p[c + 1]; j++) {
            res(groups[i[j]], c)++;
        }
    }    
    return res;
}


// [[Rcpp::export]]
std::list<float> cpp_in_place_rank_mean(arma::vec& v_temp, int idx_begin, 
                                        int idx_end) {
    std::list<float> ties;
    
    if (idx_begin > idx_end) return ties;
    std::vector<pair<float, size_t> > v_sort(idx_end - idx_begin + 1);
    for (size_t i = idx_begin; i <= idx_end; i++) {
        v_sort[i - idx_begin] = make_pair(v_temp[i], i - idx_begin);
    }
    
    
    sort(v_sort.begin(), v_sort.end());

    float rank_sum = 0, n = 1;    
    size_t i;
    for (i = 1U; i < v_sort.size(); i++) {
        if (v_sort[i].first != v_sort[i - 1].first) {
            // if current val != prev val
            // set prev val to something
            for (unsigned j = 0; j < n; j++) {
                v_temp[v_sort[i - 1 - j].second + idx_begin] = 
                    (rank_sum / n) + 1;  
            }            
            // restart count ranks
            rank_sum = i;
            if (n > 1) ties.push_back(n);
            n = 1;
        } else {
            // if curr val is a tie, 
            // don't set anything yet, start computing mean rank
            rank_sum += i;
            n++;
        }
    }
    // set the last element(s), and record the final tie group, which
    // used to be dropped from the variance correction (#29)
    for (unsigned j = 0; j < n; j++)
        v_temp[v_sort[i - 1 - j].second + idx_begin] = (rank_sum / n) + 1;
    if (n > 1) ties.push_back(n);

    return ties;
}


// [[Rcpp::export]]
std::vector<std::list<float> > cpp_rank_matrix_dgc(
        arma::vec& x, const arma::vec& p, int nrow, int ncol) {
    vector<list<float> > ties(ncol);
    int n_zero;
    for (int i = 0; i < ncol; i++) {
        n_zero = nrow - (p[i+1] - p[i]);
        if (p[i+1] == p[i]) {
            // all-zero column: the zeros are one big tie group, which used
            // to be dropped from the variance correction (#29)
            ties[i].push_back(n_zero);
            continue;
        }
        ties[i] = cpp_in_place_rank_mean(x, p[i], p[i + 1] - 1);
        ties[i].push_back(n_zero);
        x.rows(p[i], p[i + 1] - 1) += n_zero;
    }
    return ties;
}


// [[Rcpp::export]]
Rcpp::List cpp_rank_matrix_dense(const arma::mat& X_in) {
    // Rank the transposed input per column. Work on an owned copy: this
    // function used to take a non-const reference that aliased R's memory
    // and overwrote the caller's matrix with ranks (#7).
    arma::mat X = X_in.t();
    // sizes of tied groups
    vector<list<float> > ties(X.n_cols);

    std::vector<pair<float, size_t> > v_sort(X.n_rows);
    for (unsigned c = 0; c < X.n_cols; c++) {
        for (size_t i = 0; i < X.n_rows; i++) {
            v_sort[i] = make_pair(X.col(c)[i], i);
        }
        sort(v_sort.begin(), v_sort.end());

        float rank_sum = 0, n = 1;
        size_t i;
        for (i = 1U; i < v_sort.size(); i++) {
            if (v_sort[i].first != v_sort[i - 1].first) {
                // if current val != prev val
                // set prev val to something
                for (unsigned j = 0; j < n; j++) {
                    X.col(c)[v_sort[i - 1 - j].second] = (rank_sum / n) + 1;
                }
                // restart count ranks
                rank_sum = i;
                if (n > 1) ties[c].push_back(n);
                n = 1;
            } else {
                // if curr val is a tie,
                // don't set anything yet, start computing mean rank
                rank_sum += i;
                n++;
            }
        }
        // set the last element(s), and record the final tie group, which
        // used to be dropped from the variance correction (#29)
        for (unsigned j = 0; j < n; j++)
            X.col(c)[v_sort[i - 1 - j].second] = (rank_sum / n) + 1;
        if (n > 1) ties[c].push_back(n);
    }
    return Rcpp::List::create(Named("X_ranked") = X, Named("ties") = ties);
}




// [[Rcpp::export]]
arma::mat cpp_nnzeroGroups_dgc_T(const arma::vec& p, const arma::vec& i, 
                                 int ncol, int nrow, const arma::uvec& groups,
                                 int ngroups) {
    // Here, columns are samples
    arma::mat res = arma::zeros<arma::mat>(ngroups, nrow);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            // group_map gives the group num of that row
            res(groups[c], i[j])++;
        }
    }
    return res;
}


// Fused group-wise sum and non-zero count over a dense matrix in one pass.
// X is features x observations; observations are columns (MARGIN = 1).
// [[Rcpp::export]]
Rcpp::List cpp_sumGroups_nnz_dense_T(const arma::mat& X,
                                     const arma::uvec& groups,
                                     unsigned ngroups) {
    arma::mat sums = arma::zeros<arma::mat>(ngroups, X.n_rows);
    arma::mat nnz = arma::zeros<arma::mat>(ngroups, X.n_rows);
    for (unsigned c = 0; c < X.n_cols; c++) {
        unsigned g = groups[c];
        for (unsigned r = 0; r < X.n_rows; r++) {
            double v = X(r, c);
            sums(g, r) += v;
            if (v != 0) nnz(g, r)++;
        }
    }
    return Rcpp::List::create(Rcpp::Named("sums") = sums,
                              Rcpp::Named("nnz") = nnz);
}


// Compute every per-group sparse statistic wilcoxauc() needs in a single
// traversal of a dgCMatrix. X is features (nfeature) x observations (ncell) in
// CSC layout (p over observations, i = feature index, x = value). This folds
// what were five separate passes -- the R-level Matrix::t(), the per-feature
// ranking, sumGroups() over the ranked matrix, and sumGroups()/nnzeroGroups()
// over the original -- into one transpose plus one per-feature ranking.
//
// Ranking matches cpp_rank_matrix_dgc()/cpp_in_place_rank_mean(): the
// stored (non-zero) values are averaged-rank ranked amongst themselves,
// then shifted up by the number of implicit zeros in the feature. `ties`
// records every tie group per feature, including the final group among the
// stored values and the implicit-zero run, so the Wilcoxon variance
// correction is complete (#29).
// Returns `grs` (rank sums), `sums` (raw sums), `nnz` (non-zero counts) --
// each ngroups x nfeature -- and `ties` (list per feature).
//
// The per-feature ranking loop is embarrassingly parallel (each feature writes
// its own disjoint column of `grs` and its own `ties` slot), so it is split
// across `nthreads` std::threads with no locks. `nthreads = 1` runs serially.
// No R/Rcpp API is touched inside the threads: ties are accumulated into plain
// C++ vectors and converted to the R list afterwards.
// [[Rcpp::export]]
Rcpp::List cpp_wilcox_stats_dgc(const arma::vec& x, const arma::vec& p,
                                const arma::uvec& i, int nfeature, int ncell,
                                const arma::uvec& groups, int ngroups,
                                int nthreads = 1) {
    long nstored = x.n_elem;

    arma::mat sums = arma::zeros<arma::mat>(ngroups, nfeature);
    arma::mat nnz = arma::zeros<arma::mat>(ngroups, nfeature);

    // 1. Counting-sort transpose. The same pass accumulates the raw group
    //    sums and non-zero counts (used for the U statistic, pct_in/pct_out,
    //    avgExpr and logFC), so the matrix is only traversed once.
    std::vector<int> fp(nfeature + 1, 0);
    for (long k = 0; k < nstored; k++) fp[i[k] + 1]++;
    for (int f = 0; f < nfeature; f++) fp[f + 1] += fp[f];

    std::vector<unsigned> tgroup(nstored);   // group of each stored entry
    std::vector<double> tval(nstored);       // value of each stored entry
    std::vector<int> cursor(fp.begin(), fp.begin() + nfeature);
    for (int c = 0; c < ncell; c++) {
        unsigned g = groups[c];
        for (int k = (int) p[c]; k < (int) p[c + 1]; k++) {
            int f = i[k];
            double v = x[k];
            sums(g, f) += v;
            nnz(g, f)++;
            int pos = cursor[f]++;
            tgroup[pos] = g;
            tval[pos] = v;
        }
    }

    // 2 & 3. Rank each feature and accumulate group rank sums. Ties are
    // collected per feature into plain C++ vectors (no R API in threads); the
    // R list is assembled after all threads join.
    arma::mat grs = arma::zeros<arma::mat>(ngroups, nfeature);
    std::vector<std::vector<double> > ties_vec(nfeature);

    // Rank features [f_begin, f_end); writes only to disjoint columns of grs
    // and to its own ties_vec slots, so no synchronisation is needed.
    auto rank_features = [&](int f_begin, int f_end) {
        std::vector<std::pair<double, int> > v_sort;
        for (int f = f_begin; f < f_end; f++) {
            int b = fp[f];
            int m = fp[f + 1] - b;         // stored non-zeros in this feature
            std::vector<double>& feat_ties = ties_vec[f];
            if (m == 0) {
                // all-zero feature: the zeros are one big tie group (#29)
                feat_ties.push_back(ncell);
                continue;
            }
            int n_zero = ncell - m;

            v_sort.resize(m);
            for (int t = 0; t < m; t++)
                v_sort[t] = std::make_pair(tval[b + t], t);
            std::sort(v_sort.begin(), v_sort.end());

            double rank_sum = 0;
            int n = 1, t;
            for (t = 1; t < m; t++) {
                if (v_sort[t].first != v_sort[t - 1].first) {
                    double rank = (rank_sum / n) + 1 + n_zero;
                    for (int j = 0; j < n; j++)
                        grs(tgroup[b + v_sort[t - 1 - j].second], f) += rank;
                    rank_sum = t;
                    if (n > 1) feat_ties.push_back(n);
                    n = 1;
                } else {
                    rank_sum += t;
                    n++;
                }
            }
            // Assign ranks for the final tie group and record it, then the
            // implicit-zero run (#29: both used to be dropped from the
            // variance correction).
            double rank = (rank_sum / n) + 1 + n_zero;
            for (int j = 0; j < n; j++)
                grs(tgroup[b + v_sort[t - 1 - j].second], f) += rank;
            if (n > 1) feat_ties.push_back(n);
            feat_ties.push_back(n_zero);
        }
    };

    int nthr = nthreads < 1 ? 1 : nthreads;
    if (nthr > nfeature) nthr = nfeature > 0 ? nfeature : 1;
    if (nthr == 1) {
        rank_features(0, nfeature);
    } else {
        std::vector<std::thread> pool;
        int chunk = (nfeature + nthr - 1) / nthr;
        for (int t = 0; t < nthr; t++) {
            int a = t * chunk;
            int e = std::min(nfeature, a + chunk);
            if (a >= e) break;
            pool.emplace_back(rank_features, a, e);
        }
        for (std::thread& th : pool) th.join();
    }

    Rcpp::List ties(nfeature);
    for (int f = 0; f < nfeature; f++) ties[f] = ties_vec[f];

    return Rcpp::List::create(Rcpp::Named("grs") = grs,
                              Rcpp::Named("sums") = sums,
                              Rcpp::Named("nnz") = nnz,
                              Rcpp::Named("ties") = ties);
}

