#ifndef PDCGM_SMatrix_H
#define PDCGM_SMatrix_H
#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <cassert>
#include "ipm/ipx/lp_solver.h"

#define PDCGM_BIG			  1.0e+50
using Int = ipxint;

#define PDCGM_TOL_ZERO 			1.0e-10

class PDCGM_SMatrix_CW {
public:
    Int m_m{};
    Int m_n{};
    Int m_nz{};

    Int m_max_m{};
    Int m_max_n{};
    Int m_max_nz{};

    std::vector<double> m_obj{};

    std::vector<double> m_primal_lb{};
    std::vector<double> m_primal_ub{};
    std::vector<double> m_primal_lb_fxd{};
    std::vector<double> m_primal_ub_fxd{};

    std::vector<double> m_dual_lb{};
    std::vector<double> m_dual_ub{};
    std::vector<double> m_dual_lb_fxd{};
    std::vector<double> m_dual_ub_fxd{};

    std::vector<double> m_coeff{};
    std::vector<Int>    m_rwnmbs{};
    std::vector<Int>    m_clpnts{};

    std::vector<double> m_rhs{};
    std::vector<char> m_constr_type{};


    PDCGM_SMatrix_CW(){};



    //Set functions
    void PDCGM_set_rhs(const std::vector<double>& rhs) {this->m_rhs = rhs;}

    void PDCGM_set_primal_lb(const int size, const double primal_lb){this->m_primal_lb.resize(size,primal_lb);}

    void PDCGM_set_primal_lb(const std::vector<double>& primal_lb){this->m_primal_lb = primal_lb;}

    void PDCGM_set_primal_lb_fxd(const int size, const double primal_lb){this->m_primal_lb_fxd.resize(size,primal_lb);}

    void PDCGM_set_primal_lb_fxd(const std::vector<double>& primal_lb){this->m_primal_lb_fxd = primal_lb;}

    void PDCGM_set_primal_ub(const int size, const double primal_ub=INFINITY){this->m_primal_ub.resize(size,primal_ub);}

    void PDCGM_set_primal_ub(const std::vector<double>& primal_ub){this->m_primal_ub = primal_ub;}

    void PDCGM_set_primal_ub_fxd(const int size, const double primal_ub=INFINITY){this->m_primal_ub_fxd.resize(size,primal_ub);}

    void PDCGM_set_primal_ub_fxd(const std::vector<double>& primal_ub){this->m_primal_ub_fxd = primal_ub;}

    void PDCGM_set_dual_lb(const int size, const double dual_lb){this->m_dual_lb = std::vector<double>(size,dual_lb);}

    void PDCGM_set_dual_lb(const std::vector<double>& dual_lb){this->m_dual_lb = dual_lb;}

    void PDCGM_set_dual_lb_fxd(const int size, const double dual_lb){this->m_dual_lb_fxd = std::vector<double>(size,dual_lb);}

    void PDCGM_set_dual_lb_fxd(const std::vector<double>& dual_lb){this->m_dual_lb_fxd = dual_lb;}

    void PDCGM_set_dual_ub(const int size, const double dual_ub){this->m_dual_ub = std::vector<double>(size,dual_ub);}

    void PDCGM_set_dual_ub(const std::vector<double>& dual_ub){this->m_dual_ub = dual_ub;}

    void PDCGM_set_dual_ub_fxd(const int size, const double dual_ub){this->m_dual_ub_fxd = std::vector<double>(size,dual_ub);}

    void PDCGM_set_dual_ub_fxd(const std::vector<double>& dual_ub){this->m_dual_ub_fxd = dual_ub;}

    void PDCGM_set_coeff(const std::vector<double>& coeff){this->m_coeff = coeff;}

    void PDCGM_set_rwnmbs(const std::vector<Int>& rwnbms){this->m_rwnmbs = rwnbms;}

    void PDGCGM_set_clpnts(const std::vector<Int>& clpnts){this->m_clpnts = clpnts;}

    void PDCGM_set_constr_type(const std::vector<char>& constr_type) {this->m_constr_type = constr_type;}

    void PDCGM_set_constr_type(const int i, const char value) {this->m_constr_type[i] = value;}




    Int loadModel(ipx::LpSolver& lps);

    int PDCGM_set_SMatrix_CW(const int max_m,const int max_n,const int max_nz);
   
    int PDCGM_add_SMatrix_CW(const PDCGM_SMatrix_CW& M);

    int PDCGM_add_art_SMatrix_CW(const int size);

    void setStandardPrimalBounds();

    void setStandardRowInformation(const PDCGM_SMatrix_CW& M);

    int PDCGM_increase_bounds_art_SMatrix_CW(const int size,const std::vector<double>& LB, const std::vector<double>& UB);

    void PDCGM_free_SMatrix_CW();

    int PDCGM_add_dense_columns(const std::vector<double>& columns, const int n_coeff_column, const int n_columns);
    
    int PDCGM_get_m(){return m_m;};

    void update_n(const int newN);

    void update_nz(const int newN);

    void eliminate_columns(const PDCGM_SMatrix_CW& M, const std::vector<int>& indexVec,const int n_vars_art);

    void print_SMatrix_CW(const int from=0) const;

    void print_SMatrix_CW(const int from=0);

    void copy_SMatrix_CW(const PDCGM_SMatrix_CW& M, const int n_vars);
};

#endif
