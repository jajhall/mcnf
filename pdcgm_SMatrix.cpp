#include "pdcgm_SMatrix.h"


int PDCGM_SMatrix_CW::PDCGM_set_SMatrix_CW(const int max_m,const int max_n,const int max_nz)
{
    try{
    if (max_n > 0)
    {
        //std::cout << max_n << '\n';
        m_clpnts.reserve(max_n + 1);
        m_obj.reserve(max_n);
        m_primal_lb.reserve(max_n);
        m_primal_ub.reserve(max_n);
        m_primal_lb_fxd.reserve(max_n);
        m_primal_ub_fxd.reserve(max_n);
        m_clpnts.push_back(0);
    }
    else return -1;

    if (max_nz > 0)
    {
        m_coeff.reserve(max_nz);
        m_rwnmbs.reserve(max_nz);
    }
    else return -1;
    

    if(max_m>=0)
    {
        m_dual_lb.reserve(max_m);
        m_dual_ub_fxd.reserve(max_m);
        m_dual_ub.reserve(max_m);
        m_dual_ub_fxd.reserve(max_m);
        m_rhs.reserve(max_m);
        m_constr_type.reserve(max_m);
    }

    else return -1;
    }
    catch(...)
    {
        return -1;
    }

    m_max_m = max_m;
    m_max_n = max_n;
    m_max_nz = max_nz;

    m_m = max_m;
    m_n = 0;
    m_nz = 0;
    return 0;
}

/*void PDCGM_SMatrix_CW::PDCGM_free_SMatrix_CW()
{
    if(this!=NULL)
    {
        m_m =0;
        m_n=0;
        m_nz=0;
        m_max_m=0;
        m_max_n=0;
        m_max_nz=0;
        m_coeff = std::vector<double>();
        m_rwnmbs = std::vector<INTS>();
        m_clpnts = std::vector<long>();
    }
}*/

int PDCGM_SMatrix_CW::PDCGM_add_SMatrix_CW(const PDCGM_SMatrix_CW& Mnew)
{
    //Mnew.print_SMatrix_CW();

    if(Mnew.m_n + m_n > m_max_n)
	{
		update_n(m_n+(5*Mnew.m_n)+1000);
	}
	if(Mnew.m_nz + m_nz > m_max_nz)
	{
		update_nz(m_nz + 5*Mnew.m_nz + 1000);
	}

    for(int i=0; i<Mnew.m_n; ++i)
    {
        int temp = m_nz + Mnew.m_clpnts[i+1];
        m_clpnts.push_back(temp);
        m_obj.push_back(Mnew.m_obj[i]);
        m_primal_lb.push_back(Mnew.m_primal_lb[i]);
        m_primal_lb_fxd.push_back(Mnew.m_primal_lb_fxd[i]);
        m_primal_ub.push_back(Mnew.m_primal_ub[i]);
        m_primal_ub_fxd.push_back(Mnew.m_primal_ub_fxd[i]);
    }

    for(int i=0; i<Mnew.m_nz; ++i)
    {
        m_coeff.push_back(Mnew.m_coeff[i]);
        m_rwnmbs.push_back(Mnew.m_rwnmbs[i]);
    }
    m_nz+=Mnew.m_nz;
    m_n+=Mnew.m_n;
    m_clpnts[m_n] = m_nz;

    return Mnew.m_n;
}

void PDCGM_SMatrix_CW::update_n(const int newN)
{
    this->m_max_n = newN;

    this->m_clpnts.reserve(newN);
    this->m_obj.reserve(newN);
    this->m_primal_lb.reserve(newN);
    this->m_primal_ub.reserve(newN);
}

void PDCGM_SMatrix_CW::update_nz(const int newNZ)
{
    this->m_max_nz = newNZ;
    this->m_coeff.reserve(newNZ);
    this->m_rwnmbs.reserve(newNZ);
}

void  PDCGM_SMatrix_CW::print_SMatrix_CW(const int from) const
{
    /*std::cout << "m: " << m_m << " n: " << m_n <<" nz: " << m_nz << std::endl;
    std::cout << "clpnts: \n";
    for(int i=from; i<m_n+1; ++i)
    {
        std::cout << m_clpnts[i] << "  ";
    }
    std::cout << "\nCoef \trow\n";
    for(int i=m_clpnts[from]-1; i<m_nz;++i)
    {
        std::cout << m_coeff[i] <<'\t' << m_rwnmbs[i] << '\n';
    }
    std::cout << "\n\n";*/

    for(int i=0; i<m_n; ++i)
    {
        std::cout << m_obj[i] << '\t';
    }
    std::cout << '\n';
    for(int j=0; j<m_m; ++j)
    {
        for(int i=0; i<m_n; ++i)
        {
            for(int k=m_clpnts[i]; k<m_clpnts[i+1];++k)
            {
                if(m_rwnmbs[k]>j) break;
                else if(m_rwnmbs[k]==j)
                {
                    std::cout << m_coeff[k];
                }
            }
            std::cout << '\t';
        }
        std::cout << m_constr_type[j] << " " << m_rhs[j] << '\n';
    }
}

void  PDCGM_SMatrix_CW::print_SMatrix_CW(const int from)
{
    /*std::cout << "m: " << m_m << " n: " << m_n <<" nz: " << m_nz << std::endl;
    std::cout << "clpnts: \n";
    for(int i=from; i<m_n+1; ++i)
    {
        std::cout << m_clpnts[i] << "  ";
    }
    std::cout << "\nCoef \trow\n";
    for(int i=m_clpnts[from]; i<m_nz;++i)
    {
        std::cout << m_coeff[i] <<'\t' << m_rwnmbs[i] << '\n';
    }
    std::cout << "\n\n";*/
    std::cout << "Start printing matrix\n\n";
    for(int i=0; i<m_n; ++i)
    {
        std::cout << m_obj[i] << '\t';
    }
    std::cout << '\n';
    for(int j=0; j<m_m; ++j)
    {
        for(int i=0; i<m_n; ++i)
        {
            for(int k=m_clpnts[i]; k<m_clpnts[i+1];++k)
            {
                if(m_rwnmbs[k]>j) break;
                else if(m_rwnmbs[k]==j)
                {
                    std::cout << m_coeff[k];
                }
            }
            std::cout << '\t';
        }
        std::cout << m_constr_type[j] << " " << m_rhs[j] << '\n';
    }
    std::cout << "\nFinished printing matrix\n";
}

int PDCGM_SMatrix_CW::PDCGM_add_art_SMatrix_CW(const int size)
{
    PDCGM_SMatrix_CW newMatrix;
    newMatrix.PDCGM_set_SMatrix_CW(this->m_max_m,size,size*2);
    for(int i=0;i <size;++i)
    {
        newMatrix.m_coeff.push_back(1.0);
        newMatrix.m_rwnmbs.push_back(i);
        newMatrix.m_obj.push_back(this->m_dual_ub[i]);


        ++(newMatrix.m_n);
        ++(newMatrix.m_nz);

        newMatrix.m_clpnts.push_back(newMatrix.m_nz);
    }

    for(int i=0; i<size;++i)
    {   
        newMatrix.m_coeff.push_back(-1.0);
        newMatrix.m_rwnmbs.push_back(i);
        newMatrix.m_obj.push_back(-this->m_dual_lb[i]);


        ++(newMatrix.m_n);
        ++(newMatrix.m_nz);

        newMatrix.m_clpnts.push_back(newMatrix.m_nz);
    }
    newMatrix.setStandardPrimalBounds();

    PDCGM_add_SMatrix_CW(newMatrix);
    std::cout << "AFTER artificial" << '\n';
    //print_SMatrix_CW();
    return 0;
}

void PDCGM_SMatrix_CW::setStandardPrimalBounds()
{
    PDCGM_set_primal_lb(m_n,0.0);
    PDCGM_set_primal_lb_fxd(m_n,0.0);
    PDCGM_set_primal_ub(m_n,INFINITY);
    PDCGM_set_primal_ub_fxd(m_n,INFINITY);
}

void PDCGM_SMatrix_CW::setStandardRowInformation(const PDCGM_SMatrix_CW& M)
{
    this->m_rhs = M.m_rhs;
    this->m_constr_type = M.m_constr_type;
    this->m_dual_lb = M.m_dual_lb;
    this->m_dual_lb_fxd = M.m_dual_lb_fxd;
    this->m_dual_ub = M.m_dual_ub;
    this->m_dual_ub_fxd = M.m_dual_ub_fxd;
}

int PDCGM_SMatrix_CW::PDCGM_increase_bounds_art_SMatrix_CW(const int size,const std::vector<double>& LB, const std::vector<double>& UB)
{
    int aux;
    for(int i=0;i<size;++i)
    {
        aux=2*i+1;
        m_coeff[aux] = UB[i];

        aux=2*i+2*size+1;
        m_coeff[aux] = -LB[i];
    }
    return 0;
}

int PDCGM_SMatrix_CW::PDCGM_add_dense_columns(const std::vector<double>& columns, const int n_coeff_column, const int n_columns)
   {
        PDCGM_SMatrix_CW newMatrix;
        newMatrix.PDCGM_set_SMatrix_CW(this->m_max_m,n_columns,n_columns*n_coeff_column);

        for(int i=0; i< n_columns;++i)
        {
            newMatrix.m_obj.push_back(columns[n_coeff_column*i]);
            for(int j=1; j<n_coeff_column;++j)
            {
                double tempVal = columns[n_coeff_column*i+j];
                if(fabs(tempVal)>PDCGM_TOL_ZERO)
                {
                    newMatrix.m_coeff.push_back(tempVal);
                    newMatrix.m_rwnmbs.push_back(j-1);
                    ++(newMatrix.m_nz);
                }
            }
            ++(newMatrix.m_n);
            newMatrix.m_clpnts.push_back(newMatrix.m_nz);
        }
        if(newMatrix.m_nz==0) return 0;
        newMatrix.setStandardPrimalBounds();
        PDCGM_add_SMatrix_CW(newMatrix);
        std::cout << "After adding dense column\n";
        //print_SMatrix_CW();
        return newMatrix.m_n;
   } 
    

void PDCGM_SMatrix_CW::eliminate_columns(const PDCGM_SMatrix_CW& M, const std::vector<int>& indexVec, const int n_vars_art)
{
    
    //assert(1==0);
    int column;
    assert(m_n==0);
    std::cout << "Eliminating columns" << std::endl;
    printf("n: %ld m: %ld nz:%ld\n",M.m_n,M.m_m,M.m_nz);
    printf("size of old rwnmbs %ld \t size of old coeff:%ld \t old nz: %ld \t",M.m_rwnmbs.size(),M.m_coeff.size(),M.m_nz);
    printf("last column thing: %ld",M.m_clpnts[M.m_n]);
    std::cout << std::endl;
    for(int i=0; i<n_vars_art;++i)
    {
        m_nz+= M.m_clpnts[m_n+1]-M.m_clpnts[m_n];
        m_clpnts.push_back(m_nz);
        m_obj.push_back(M.m_obj[m_n]);
        for(int j=M.m_clpnts[m_n]; j<M.m_clpnts[m_n+1]; ++j)
        {
            m_rwnmbs.push_back(M.m_rwnmbs[j]);
            m_coeff.push_back(M.m_coeff[j]);
        }
        ++m_n;
    }
    for(int i=0; i< indexVec.size(); ++i)
    {
        if(indexVec[i]==-1) break;
        else
        {
            column = indexVec[i] + n_vars_art;
            ++(m_n);
            m_nz+= M.m_clpnts[column+1]-M.m_clpnts[column];
            m_clpnts.push_back(m_nz);
            m_obj.push_back(M.m_obj[column]);
            for(int j=M.m_clpnts[column]; j<M.m_clpnts[column+1]; ++j)
            {
                m_rwnmbs.push_back(M.m_rwnmbs[j]);
                m_coeff.push_back(M.m_coeff[j]);
            }
        }
    }
    setStandardPrimalBounds();
    setStandardRowInformation(M);
}



void PDCGM_SMatrix_CW::copy_SMatrix_CW(const PDCGM_SMatrix_CW& M, const int n_vars)
{
    printf("COPYING\n");
    //M.print_SMatrix_CW(0);
    printf("\n\n\n");
    this->m_m = M.m_m;
    this->m_n = M.m_n-n_vars;
    this->m_nz = M.m_nz-n_vars;
    this->m_max_m = M.m_max_m;
    this->m_max_n = M.m_max_n;
    this->m_max_nz= M.m_max_nz;

    for(int i=n_vars; i<=M.m_n;++i) 
    {
        this->m_clpnts.push_back(M.m_clpnts[i]);
        this->m_primal_lb.push_back(M.m_primal_lb[i]);
        this->m_primal_lb_fxd.push_back(M.m_primal_lb_fxd[i]);
        this->m_primal_ub.push_back(M.m_primal_ub[i]);
        this->m_primal_ub_fxd.push_back(M.m_primal_ub_fxd[i]);
        this->m_obj.push_back(M.m_obj[i]);
    }
    
    this->m_dual_lb = M.m_dual_lb;
    this->m_dual_lb_fxd = M.m_dual_lb_fxd;
    this->m_dual_ub = M.m_dual_ub;
    this->m_dual_ub_fxd = M.m_dual_ub_fxd;
    this->m_rhs = M.m_rhs;
    this->m_constr_type = M.m_constr_type;
    this->m_rwnmbs = M.m_rwnmbs;
    this->m_coeff = M.m_coeff;
    printf("\n\n\n");
}

Int PDCGM_SMatrix_CW::loadModel(ipx::LpSolver& lps)
{
    //print_SMatrix_CW();
    /*std::cout << "Number of variables \t"<< m_n << "\nobjective \t\t";
    for(int i=0; i<m_obj.size();++i) std::cout << m_obj.data()[i] << ' ';
    std::cout << "\nprimal lb \t\t";

    for(int i=0; i<m_primal_lb.size();++i) std::cout << m_primal_lb.data()[i] << ' ';
    std::cout << "\nprimal ub \t\t";

    for(int i=0; i<m_primal_ub.size();++i) std::cout << m_primal_ub.data()[i] << ' ';
    std::cout << "\n number constraints \t" << m_m << "\nclpnts \t\t\t";

    for(int i=0; i<m_clpnts.size();++i) std::cout << m_clpnts.data()[i] << ' ';
    std::cout << "\nrow numbers\t\t";

    for(int i=0; i<m_rwnmbs.size(); ++i) std::cout << m_rwnmbs.data()[i] << ' ';
    std::cout << "\ncoeff\t\t\t";

    for(int i=0; i<m_coeff.size();++i) std::cout << m_coeff.data()[i] << ' ';
    std::cout << "\nrhs\t\t\t";

    for(int i=0; i<m_rhs.size();++i) std::cout << m_rhs.data()[i] << ' ';
    std::cout << "\nConstraint type \t";

    for(int i=0; i<m_constr_type.size(); ++i) std::cout <<m_constr_type.data()[i] << ' ';
    std::cout << "\n";*/


    //Test to check whether model is feasible

    //std::cout << "m_n : " << m_n << " m_m: " << m_m << '\n';
    //Objective is feasible
    assert(m_obj.size()==m_n);
    for(int i=0; i<m_obj.size();++i)
    {
        //std::cout << m_obj[i] << std::endl;
        if(m_obj[i]<0) assert(0);
        //assert(m_obj[i]>=0);
    }
    //std::cout << '\n';

    //Primal_LB is feasible
    assert(m_primal_lb.size()==m_n);
    for(int i=0; i<m_primal_lb.size(); ++i) 
    {
        //std::cout << m_primal_lb[i] << ' '; 
        //assert(m_primal_lb[i]==0);
    }
    //std::cout << '\n';
    

    //Primal_UB is feasible
    assert(m_primal_ub.size()==m_n);
    for(int i=0; i<m_primal_ub.size(); ++i) 
    {
        //std::cout << m_primal_ub[i] << ' '; 
        assert(m_primal_ub[i]==INFINITY);
    }
    //std::cout << '\n';


    //Clpnts is feasible
    assert(m_clpnts.size()==m_n+1);
    for(int i=0; i<m_clpnts.size();++i)
    {
        //std::cout << m_clpnts[i] << ' ';
        assert(m_clpnts[i]>=0);
        if(i>0) assert(m_clpnts[i]>=m_clpnts[i-1]);
    }
    //std::cout << '\n';

    //RWNMBS is feasible
    assert(m_rwnmbs.size()==m_nz);
    //std::cout << m_rwnmbs.size() << '\n';
    //std::cout << "m_m: " << m_m <<'\n';
    for(int i=0; i<m_rwnmbs.size(); ++i)
    {
        //std::cout <<  m_rwnmbs[i] << ' ';
        //if(m_rwnmbs[i]>=m_m) break;
        assert(m_rwnmbs[i]>=0 && m_rwnmbs[i]<m_m && i+1);
    }
   // std::cout << std::endl;
    //assert(0);
    //assert(0);

    //coeff is feasible
    assert(m_coeff.size()==m_nz);
    for(int i=0; i<m_coeff.size(); ++i)
    {
        //std::cout << m_coeff[i] << ' ';
    }
   // std::cout << std::endl;

    //std::cout << m_rhs.size() << ' ' << m_m << std::endl;
    assert(m_rhs.size()==m_m);
    for(int i=0; i<m_rhs.size(); ++i)
    {
        //std::cout << m_rhs[i] << ' ';
        
    }
    //std::cout << std::endl;

    assert(m_constr_type.size()==m_m);
    for(int i=0; i<m_constr_type.size(); ++i)
    {
        //std::cout << m_constr_type[i] << ' ';
    }

    return lps.LoadModel(m_n,m_obj.data(), m_primal_lb.data(),m_primal_ub.data(),m_m,m_clpnts.data(),m_rwnmbs.data(),m_coeff.data(),m_rhs.data(),m_constr_type.data());
}
    