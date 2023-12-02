// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/exdex.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// arma_rowSums_minus_col
arma::colvec arma_rowSums_minus_col(const arma::mat& x, const int& j);
static SEXP _exdex_arma_rowSums_minus_col_try(SEXP xSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int& >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(arma_rowSums_minus_col(x, j));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _exdex_arma_rowSums_minus_col(SEXP xSEXP, SEXP jSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_exdex_arma_rowSums_minus_col_try(xSEXP, jSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// cpp_log0const
arma::mat cpp_log0const(const arma::mat& x, const double& constant);
static SEXP _exdex_cpp_log0const_try(SEXP xSEXP, SEXP constantSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type constant(constantSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_log0const(x, constant));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _exdex_cpp_log0const(SEXP xSEXP, SEXP constantSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_exdex_cpp_log0const_try(xSEXP, constantSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// cpp_col_ms
arma::rowvec cpp_col_ms(arma::mat const& x);
static SEXP _exdex_cpp_col_ms_try(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_col_ms(x));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _exdex_cpp_col_ms(SEXP xSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_exdex_cpp_col_ms_try(xSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// cpp_sigma2hat_dj
Rcpp::List cpp_sigma2hat_dj(const Rcpp::List& all_max, const int& b, const int& kn, const int& m, const String& bias_adjust, const String& which_dj);
static SEXP _exdex_cpp_sigma2hat_dj_try(SEXP all_maxSEXP, SEXP bSEXP, SEXP knSEXP, SEXP mSEXP, SEXP bias_adjustSEXP, SEXP which_djSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type all_max(all_maxSEXP);
    Rcpp::traits::input_parameter< const int& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int& >::type kn(knSEXP);
    Rcpp::traits::input_parameter< const int& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const String& >::type bias_adjust(bias_adjustSEXP);
    Rcpp::traits::input_parameter< const String& >::type which_dj(which_djSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_sigma2hat_dj(all_max, b, kn, m, bias_adjust, which_dj));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _exdex_cpp_sigma2hat_dj(SEXP all_maxSEXP, SEXP bSEXP, SEXP knSEXP, SEXP mSEXP, SEXP bias_adjustSEXP, SEXP which_djSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_exdex_cpp_sigma2hat_dj_try(all_maxSEXP, bSEXP, knSEXP, mSEXP, bias_adjustSEXP, which_djSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _exdex_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("arma::colvec(*arma_rowSums_minus_col)(const arma::mat&,const int&)");
        signatures.insert("arma::mat(*cpp_log0const)(const arma::mat&,const double&)");
        signatures.insert("arma::rowvec(*cpp_col_ms)(arma::mat const&)");
        signatures.insert("Rcpp::List(*cpp_sigma2hat_dj)(const Rcpp::List&,const int&,const int&,const int&,const String&,const String&)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _exdex_RcppExport_registerCCallable() { 
    R_RegisterCCallable("exdex", "_exdex_arma_rowSums_minus_col", (DL_FUNC)_exdex_arma_rowSums_minus_col_try);
    R_RegisterCCallable("exdex", "_exdex_cpp_log0const", (DL_FUNC)_exdex_cpp_log0const_try);
    R_RegisterCCallable("exdex", "_exdex_cpp_col_ms", (DL_FUNC)_exdex_cpp_col_ms_try);
    R_RegisterCCallable("exdex", "_exdex_cpp_sigma2hat_dj", (DL_FUNC)_exdex_cpp_sigma2hat_dj_try);
    R_RegisterCCallable("exdex", "_exdex_RcppExport_validate", (DL_FUNC)_exdex_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_exdex_arma_rowSums_minus_col", (DL_FUNC) &_exdex_arma_rowSums_minus_col, 2},
    {"_exdex_cpp_log0const", (DL_FUNC) &_exdex_cpp_log0const, 2},
    {"_exdex_cpp_col_ms", (DL_FUNC) &_exdex_cpp_col_ms, 1},
    {"_exdex_cpp_sigma2hat_dj", (DL_FUNC) &_exdex_cpp_sigma2hat_dj, 6},
    {"_exdex_RcppExport_registerCCallable", (DL_FUNC) &_exdex_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_exdex(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
