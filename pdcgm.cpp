#include "Highs.h"
#include "McnfOracle.h"
//#include "pdcgm_env.h"

const bool dev_run = true;

static void userColOracle(const int col_oracle_type, const char* message,
                             const HighsColOracleDataOut* data_out,
                             HighsColOracleDataIn* data_in,
                             void* user_col_oracle_data);
/*
HighsInt MCNF_linear_solve(McnfGraph * G, McnfDemand * dem, HighsInt nb_dem, HighsInt aggregated);

static short MCNF_linear_disagg_oracle(double *primal_violation,   // total violation of the generated constraints (not used in this application) 
				       double *dual_violation,     // total violation of the generated columns (the value f(y)) 
				       PDCGM_env  &MCNF,          // PDCGM environment 
				       void   *instance_data);      // instance data (not used in this example) 
*/
int main(int argc, char *argv[]) {
  // Auxiliary variables
  HighsInt nb_dem, aggregated, oriented, begin_zero, with_cost;

  // Graph representing the network
  McnfGraph G; 
  // Demand of each commodity
  McnfDemand *dem; 

  // Instance file name
  char DefltFlnme[1000];

  McnfOracle oracle;
  // Parse the command line receive as an input
  oracle.parseInput(argc, argv, &aggregated, &oriented, DefltFlnme, &begin_zero, &with_cost);

  // Summarize the input data
  printf("\nProblem filename: %s\n", DefltFlnme);

  FILE *out_file = fopen("output-mcnf-pdcgm.txt", "a");
  fprintf(out_file, "\n%s ", DefltFlnme);
  fflush(out_file); 
  fclose(out_file);

  // Read the instance data 
  oracle.read_instance(DefltFlnme, begin_zero, aggregated);

  // Set initial parameters and initialize the oracle data structure 
  oracle.set_network_Dijkstra(&G, &dem, &nb_dem, oriented);
  oracle.init_Dijkstra(&G, dem, nb_dem);
  // Solve the problem by column generation
  // MCNF_linear_solve(&G, dem, nb_dem, aggregated);


  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setColOracle(userColOracle);
  return 0;
}

static void userColOracle(const int col_oracle_type, const char* message,
                             const HighsColOracleDataOut* data_out,
                             HighsColOracleDataIn* data_in,
                             void* user_col_oracle_data) {
}

/*
HighsInt MCNF_linear_solve(McnfGraph * G, McnfDemand * dem, HighsInt nb_dem, HighsInt aggregated) {
    HighsInt dim = G->m, // number of linking constraints in the problem: number of arcs 
        nb_sub_pb,   // number of subproblems (= number of convexity constraints) 
        max_n_cols = 10000,
        max_outer = 5000, 
        k;
    
    PDCGM_env MCNF; // PDCGM environment 
    
    FILE *out_file; // output file pointer  
    
    // Set the number of subproblems according to the formulation 
    nb_sub_pb = nb_dem; // number of commodities 
    
    std::vector<double> b, lo_box, up_box;
    std::vector<char> row_type;
    double* y0;

    // Allocate memory for the MP formulation 
    b.reserve(dim + nb_sub_pb);
    row_type.reserve(dim+nb_sub_pb);
    y0 = PDCGM_ALLOC(double, dim + nb_sub_pb + 1);
    lo_box.reserve(dim);
    up_box.reserve(dim);

    //double *b  = PDCGM_ALLOC(double, dim + nb_sub_pb + 1); // RHS of the master problem 
    double *y0 = PDCGM_ALLOC(double, dim + nb_sub_pb + 1); // initial guess for the dual solution 
    double *lo_box = PDCGM_ALLOC(double, dim); // lower bounds of y 
    double *up_box = PDCGM_ALLOC(double, dim); // upper bounds of y 
    short  *row_type = PDCGM_ALLOC(short, dim + nb_sub_pb + 1); // type of each constraint in the master problem 


    // Set the linking constraints and the associated duals 
    for(k = 0; k < dim; k++) 
    {
        b.push_back(-Network->K[k]);
        y0[k] = 0.0;
        row_type.push_back('=');
        lo_box.push_back(0.0);
        up_box.push_back(10.0);
    }
    
    // Set the convexity constraints and the associated duals 
    for (k = dim; k < dim + nb_sub_pb; k++) 
    {
        b.push_back(1.0);
        y0[k] = 100.0;
        row_type.push_back('='); // REMARK: WE RESET AS EQUALITY INSIDE THE ORACLE
    }
    //std::cout << "B size " << b.size() << '\n';

    // POPULATE THE INTERNAL DATA STRUCTURE OF PDCGM  
    MCNF.PDCGM_set_data( 
        dim,  // number of kept constraints in the RMP 
        nb_sub_pb,  // number of convexity constraints in the RMP 
        (dim + nb_sub_pb + 1),  // maximum number of nonzeros in a column 
        nb_sub_pb,     // maximum number of columns generated by the oracle 
        max_outer,     // maximum number of outer iterations 
        y0, // initial guess of the dual solution (may be NULL) 
        b,  // RHS of each constraint in the RMP 
        row_type, // type of each constraint (row) in the RMP 
        lo_box,   // lower bound vector of the DUAL variables in the RMP 
        up_box,   // upper bound vector of the DUAL variables in the RMP 
        NULL);   // instance data 
                
    // Check if the PDCGM environment was created successfully 
    if(PDCGM_env == NULL)   
    {
        printf("-- ERROR: PDCGM_set_data() returned NULL --\n");
        fflush(stdout);
        return 0;
    }
            
    // Set the optimality tolerance for the column generation algorithm 
    MCNF.PDCGM_set_delta( 1.E-5);   
    
    // Set the verbose mode (how much information is printed) 
    MCNF.PDCGM_set_verbosity( 1);
    
    // Set the maximum number of columns as a large value 
    if(max_n_cols < 2 * nb_sub_pb) max_n_cols = 2 * nb_sub_pb;
    
    MCNF.PDCGM_set_max_number_columns(max_n_cols);
    //MCNF.PDCGM_set_start_from_reduce_matrix(0);
    
    // Tell PDCGM that the first iterations may be difficult 
    MCNF.PDCGM_set_use_diff_iters(1);
    
    // Set maximum CPU time 
    MCNF.PDCGM_set_max_cputime(50000.0);
    
    // Enable column elimination on *insertion* 
    MCNF.PDCGM_set_column_elimination(1);
    
    
    // Set an (artificial) initial column to avoid infeasibility 
    //double *initial_column = PDCGM_ALLOC(double, dim + nb_sub_pb + 1);
    //for(k=dim; k < dim + nb_sub_pb; k++) initial_column[k] = 1.0;
    //initial_column[dim + nb_sub_pb] = 1.E+2;    
    //PDCGM_add_dense_columns(PDCGM_env, initial_column, NULL, dim + nb_sub_pb + 1, 1);

    // START THE COLUMN GENERATION PROCEDURE 
    
    // Solve the LINEAR multicommodity network flow problem 
    // Set the degree of optimality (parameter D > 1) 
    MCNF.PDCGM_set_degree_of_optimality(10.0);

    // Call PDCGM 
    MCNF.PDCGM_solve_MP( MCNF_linear_disagg_oracle);
    
    // Print the solution 
    //std::vector<double> lambda;
    //MCNF.PDCGM_get_pointer_to_master_solution(lambda);
    //if(lambda.empty()) std::cout << "ERROR obtaining primal solution variables\n";
    
    printf("Active sets: %4.2lf \n", (1.0 * active_set_n) / G->m);
    printf("Optimal value: %1.5E \n\n", max_cap*max_cost*MCNF.PDCGM_get_obj_UB());
    
    // Printf the results in an output file 
    //out_file = fopen("output-mcnf-pdcgm.txt", "a");
    fprintf(out_file, "\t%1.6E \t%1.6E \t%1.6E \t%lf \t%d \t%lf \t%lf \t%lf ",  
        max_cap*max_cost*PDCGM_env->lowerBnd, max_cap*max_cost*PDCGM_env->upperBnd, 
        PDCGM_env->rel_gap, (100.0 * active_set_n) / G->m, PDCGM_env->outer, PDCGM_env->cputime_RMP, PDCGM_env->cputime_oracle, PDCGM_env->cputime_CG);
    fflush(out_file);
    fclose(out_file);
    
    // Clean the memory before returning 
    
    return 0;

}

static short MCNF_linear_disagg_oracle(
    double *primal_violation,   // total violation of the generated constraints (not used in this application) 
    double *dual_violation,     // total violation of the generated columns (the value f(y)) 
    PDCGM_env  &MCNF,          // PDCGM environment 
    void   *instance_data)      // instance data (not used in this example) 
{
    
    HighsInt i, j, nr, nc, nnz, ret;
        
    double aux_double, cost, answer_type, *d; // a pointer to the dual solution of the current master problem 
        
    PDCGM_SMatrix_CW M{}, M2{}; // Auxiliary matrix to store the generated columns (sparse representation) 
    //PriorityQueue_type *priority_queue; // priority queue used to sort elements 
    
    printf ("Hello from MCNF_linear_disagg_oracle() -- Linear multicommodity network flow with disaggregated MP.\n");
    //std::cout << "OUTER ITERATIONS: " <<MCNF.PDCGM_get_outer() << std::endl;
    
    // Set the answer type: column generation 
    if(MCNF.PDCGM_get_outer()> 0) answer_type = 1; 
    else answer_type = 0; 
    
    // Reset the number of columns added in the last interation 
    MCNF.PDCGM_set_ncols_last(0);
    
    // Reset the pricing subproblem solver 
    Reset_Dijkstra (Network);
    
    // Set the number of dense entries in each generated column 
    nr = Network->m + nb_demande; // number of rows = number of arcs + convexity rows + objective row 
    
    // Set the max number of columns to be generated 
    nc = nb_demande;
    
    // Set the max number of nonzeros in a column 
    nnz = Network->m + 1;
    
    // Reset the dual_violation (the value of the oracle) 
    *dual_violation = 0.0;
    
    // Get a pointer to the dual solution of the current master problem 
    std::vector<double> u = MCNF.PDCGM_get_dual();

    //std::cout << " size of u: " << u.size();
    //for(HighsInt i=0; i<u.size();++i) std::cout << u[i] << '\n';
    //std::cout << "Finish printing u\n";
    
    
    // Set the network costs using the dual solution and the linear coefficient costs  
    for (j = 0; j < Network->m; j++) 
    { 
        net_cost[j] = u[j] + Network->cost[j];
        //std::cout << net_cost[j] << '\n';
    }
    //std::cout << "net cost" << '\n';

    // Set the sparse matrix that will store the generated columns 
    M.PDCGM_set_SMatrix_CW(nr, nc, nc*nnz);
    M.m_m = nr; // number of rows 
    M.m_n = 0;  // no columns yet 
    M.m_nz = 0; // no coefficient yet  
    
    
    // Set the sparse matrix that will store the generated columns 
    M2.PDCGM_set_SMatrix_CW(nr, nc, nc*nnz);
    
    // Allocate a memory for a queue that will be used to sort the values 
    std::priority_queue<std::pair<HighsInt,HighsInt>> pq{};
    
    // Solve the pricing subproblems of the first type (one shortest path for each commodity) 
    #pragma omp parallel for private(i)
    for (i = 0; i < nb_demande; i++)
    {        
        // Call the Dijkstra algorithm to solve the pricing subproblem 
        if(dijkstra_from_node[(demande + i)->s] == 0) Dijkstra_solve(Network, (demande + i)->s, (demande + i)->t, net_cost);
    }
    
    //std::cout << "finished solving the pricing problem\n";
    // Solve the pricing subproblems of the first type (one shortest path for each commodity) 
    for (i = 0; i < nb_demande; i++)
    {        
        // Set the vector of solution costs 
        d = d_tableau + (demande + i)->s * Network->n;
        
        //std::cout << d[(demande + i)->t] << ' ' <<  (demande + i)->val << ' ' << u[Network->m + i] << '\n';
        if( MCNF.PDCGM_get_outer()<1 || d[(demande + i)->t] * (demande + i)->val - u[Network->m + i] < -MCNF_TOL_ZERO)
        {
            // Generate the column from the optimal path 
            //printf("(demnde + %d)->val: %lf\n",i,(demande+i)->val);
            get_solution_from_Dijkstra(Network, &M, (demande + i)->s, (demande + i)->t, (demande + i)->val);
            //printf("After Dijkstra, nz: %d\n",M.m_nz);
            //M.print_SMatrix_CW();
            // Set the entry in the convexity row 
            M.m_coeff.push_back(1.0);
            short val = Network->m+i;
            M.m_rwnmbs.push_back( val);
            M.m_nz++;
           // printf("Setting first entry\n");
            //.print_SMatrix_CW();
            
            // Set the objective function cost 
            cost = 0.0;
            //std::cout << "Printing coeff\n";
            //for(auto i: M.m_coeff) std::cout << i << ' ';
            //std::cout << "Finsihed printing coeff "<< std::endl;
            for (j = M.m_clpnts[M.m_n]; j < M.m_coeff.size()-1; j++)
            {
                //std::cout << "j: " << j << "\tM.m_rwnmbs[j] -1 =" << M.m_rwnmbs[j] <<"\tNetwork cost " <<
                            //Network->cost[M.m_rwnmbs[j]] << "\tCoeff: "<< M.m_coeff[j] << std::endl;
                cost += M.m_coeff[j] * Network->cost[M.m_rwnmbs[j]];
                //std::cout << "cost: " <<cost << std::endl;
                //assert(M.m_coeff[j]<=0);
                //assert(Network->cost[M.m_rwnmbs[j]-1]<=0);
            }
            assert(cost<0);
            //assert(i<10);
            M.m_obj.push_back(-cost);
          //  printf("Setting second entry\n");
           // M.print_SMatrix_CW();

            pq.push(std::make_pair(d[(demande + i)->t] * (demande + i)->val - u[Network->m + i], M.m_n));
            //PriorityQueue_insert(priority_queue, d[(demande + i)->t] * (demande + i)->val - u[Network->m + i], M.m_n);

            // Update the number of columns and set the first index of the next column 
            M.m_n++;
            M.m_clpnts.push_back(M.m_nz);
        }
        
        // Set the relative cost (the value of the oracle)  
        //printf("condition LHS: %lf, RHS: %lf, \n",d[(demande + i)->t] * (demande + i)->val, u[Network->m + i]);
        if(d[(demande + i)->t] * (demande + i)->val - u[Network->m + i] < 0) 
        {
            //f("added to dual violation %lf",d[(demande + i)->t] * (demande + i)->val - u[Network->m + i]);
            *dual_violation += d[(demande + i)->t] * (demande + i)->val - u[Network->m + i];
        }
        if(MCNF.PDCGM_get_rwstat(Network->m+i) == '>') MCNF.PDCGM_set_rwstat(Network->m+i,'=');
    }
    //assert(0);
    //M.print_SMatrix_CW();
   // std::cout << '\n';
    while( ( (MCNF.PDCGM_get_outer() < 1) || (M2.m_n < nb_demande) ) && (!pq.empty()) )
    {
        std::pair<HighsInt,HighsInt> aux = pq.top();
        pq.pop();
        i = aux.second;
        for (j = M.m_clpnts[i]; j < M.m_clpnts[i+1]; j++)
        {
            //printf("j: %d\n",j);
            M2.m_coeff.push_back(M.m_coeff[j]);
            M2.m_rwnmbs.push_back(M.m_rwnmbs[j]);
            M2.m_nz++;
        }
        M2.m_obj.push_back(M.m_obj[i]);
        M2.m_n++;
        M2.m_clpnts.push_back(M2.m_nz);
    }

        
    //M2.print_SMatrix_CW();
    // Add the new columns to the current master problem 
    //M2.print_SMatrix_CW();
    //assert(0);
    M2.setStandardPrimalBounds();
    nc = MCNF.PDCGM_add_columns(M2, NULL);
    //printf("\n\n\n\n");
    //MCNF.PDCGM_print_Matrix(0);

    
    //assert(nc==0);
    //printf("before leaving oracle, nc=%d and dual_violation=%lf\n",nc,*dual_violation);
    if(nc == 0) *dual_violation = 0.0;
    
    return answer_type;
}
*/
