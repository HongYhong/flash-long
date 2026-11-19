#include <RcppArmadillo.h>

using namespace Rcpp;

static List em_theta_impl(const arma::mat &X,
                          const arma::rowvec &Y,
                          const int maxiter,
                          const double minvalue,
                          const double conv) {
    int M = static_cast<int>(X.n_rows);
    arma::rowvec theta(M);
    theta.fill(1.0);

    arma::mat theta_trace(M, maxiter);
    int iter = 0;
    theta_trace.col(iter) = theta.t();

    arma::mat Xb = X;
    arma::mat lmat = X;
    arma::rowvec summed_by_row_Xb = arma::sum(Xb.t(), 0);
    arma::mat adjP_Xt = X;
    arma::rowvec summed_by_col_adjP_Xt = arma::sum(adjP_Xt, 0);

    double deltaTheta = 1.0;
    arma::vec t_before = theta_trace.col(iter);
    arma::vec t_after = theta_trace.col(iter);
    arma::vec deltaVec(M);
    while (deltaTheta > conv && iter < (maxiter - 1)) {
        iter++;
        adjP_Xt = (X.t() * diagmat(theta)).t();
        summed_by_col_adjP_Xt = arma::sum(adjP_Xt, 0);
        lmat = (adjP_Xt * diagmat(Y / summed_by_col_adjP_Xt));
        lmat.replace(arma::datum::nan, 0.0);
        theta = arma::sum(lmat.t(), 0) / summed_by_row_Xb;
        theta_trace.col(iter) = theta.t();
        t_before = theta_trace.col(iter - 1);
        t_after = theta_trace.col(iter);
        deltaVec = abs(t_after.elem(find(t_after > minvalue)) -
                       t_before.elem(find(t_after > minvalue))) /
                   t_after.elem(find(t_after > minvalue));
        if (deltaVec.is_empty()) {
            deltaTheta = 0.0;
        } else {
            deltaTheta = max(deltaVec);
        }
    }

    List ret;
    ret["theta"] = theta;
    ret["theta_trace"] = theta_trace;
    ret["iter"] = iter;
    return ret;
}

extern "C" SEXP bambu_next_em_theta(SEXP XSEXP, SEXP YSEXP, SEXP maxiterSEXP, SEXP minvalueSEXP, SEXP convSEXP) {
    BEGIN_RCPP
    const arma::mat X = Rcpp::as<arma::mat>(XSEXP);
    const arma::rowvec Y = Rcpp::as<arma::rowvec>(YSEXP);
    const int maxiter = Rcpp::as<int>(maxiterSEXP);
    const double minvalue = Rcpp::as<double>(minvalueSEXP);
    const double conv = Rcpp::as<double>(convSEXP);
    return Rcpp::wrap(em_theta_impl(X, Y, maxiter, minvalue, conv));
    END_RCPP
}

extern "C" SEXP bambu_next_emWithL1(SEXP ASEXP, SEXP A_fullSEXP, SEXP A_uniqueSEXP,
                                     SEXP YSEXP, SEXP KSEXP,
                                     SEXP maxiterSEXP, SEXP minvalueSEXP, SEXP convSEXP) {
    BEGIN_RCPP
    const arma::mat A = Rcpp::as<arma::mat>(ASEXP);
    const arma::mat A_full = Rcpp::as<arma::mat>(A_fullSEXP);
    const arma::mat A_unique = Rcpp::as<arma::mat>(A_uniqueSEXP);
    const arma::rowvec Y = Rcpp::as<arma::rowvec>(YSEXP);
    const arma::rowvec K = Rcpp::as<arma::rowvec>(KSEXP);
    const int maxiter = Rcpp::as<int>(maxiterSEXP);
    const double minvalue = Rcpp::as<double>(minvalueSEXP);
    const double conv = Rcpp::as<double>(convSEXP);

    int M = static_cast<int>(A.n_rows);
    List theta_out = em_theta_impl(A, Y, maxiter, minvalue, conv);
    arma::rowvec theta = Rcpp::as<arma::rowvec>(theta_out[0]);

    arma::mat estMat(3, M);
    arma::rowvec baseSum = K % Y / arma::sum((A.t() * diagmat(theta)).t(), 0);
    baseSum.replace(arma::datum::nan, 0.0);
    estMat.row(0) = arma::sum(((A.t() * diagmat(theta)).t() * diagmat(baseSum)).t(), 0);
    estMat.row(1) = arma::sum(((A_full.t() * diagmat(theta)).t() * diagmat(baseSum)).t(), 0);
    estMat.row(2) = arma::sum(((A_unique.t() * diagmat(theta)).t() * diagmat(baseSum)).t(), 0);

    List ret;
    ret["theta"] = estMat;
    return Rcpp::wrap(ret);
    END_RCPP
}
