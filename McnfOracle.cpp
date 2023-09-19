/*
 *  PDCGM  (Primal-Dual Column Generation Method)
 *  Copyright:     Jacek Gondzio, 1996, 2013.
 *
 *  COLUMN GENERATION FOR MULTICOMMODITY NETWORK FLOW PROBLEMS
 *  USING TWO DIFFERRENT FORMULATIONS: AGGREGATED AND DISAGGREGATED
 * 
 *  Authors: Pedro Munari [munari@dep.ufscar.br], 
 *           Pablo González-Brevis, 
 *           Jacek Gondzio
 *  
 *  AUXILIARY FUNCTIONS FOR THE ORACLE PROCEDURE
 *  Last Modified: July, 2013.
 * 
 *  On top of HOPDM  (Higher Order Primal-Dual Method)
 *  Copyright:     Jacek Gondzio, 1990, 2010.
 *
 *  Using the Dijkstra algorithm implemented by:
 *      Jacek Gondzio and Robert Sarkissian
 *
 */

/* Bring in auxiliary headers */
#include "McnfOracle.h"

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <climits>// For INT_MIN
#include <cfloat>// For DBL_MAX

/* Parse the input given in the text line */
void McnfOracle::parseInput(int argc, char *argv[], HighsInt *aggregated, HighsInt *oriented, 
			    char DefltFlnme[], HighsInt *begin_zero, HighsInt *with_cost) {
  /* Default values of control parameters */
  *aggregated = 0;
  *oriented = 1;
  *begin_zero = 0;
  *with_cost = 1;
  USE_ACTIVE_SET_ = 0;
  if (argc > 1) {
      /* Get the name of the instance file */
      strcpy (DefltFlnme, argv[1]);
  } else {
    printf("\nError: Instance file name not specified.");
    std::exit(1);
  }
    
  if(argc > 2)    {
    /* Check if aggregation is required */
    if(argv[2][0] == '1') {
      *aggregated = 1;
    } else if(argv[2][0] == '2') {
      *aggregated = 2;
    }
  } 
  
  if(argc > 3)    {
    /* Check if aggregation is required */
    if(argv[3][0] == '0') 
      {
	USE_ACTIVE_SET_ = 0;
      }
  } 
}

/* Read the instance data from an input file */
void McnfOracle::read_instance(char file_name[], HighsInt zero, HighsInt aggregated) {   
  FILE *input_file;
  HighsInt    i, p;
  char   *ret, firstline[1000];// = "#Instance: \n";
  double flot, scale_factor, aux_double;

  /* Open the instance file */
  if (!(input_file = fopen(file_name, "r+t")))     {
    printf("\nError: Could not open input file %s!\n", file_name);
    exit(1);
  }

  /* Ignore the first line of the file */
  ret = fgets(firstline, 1000, input_file);
    
  /* Read the dimensions of the instance */
  ret = fgets(firstline, 1000, input_file);
  p = sscanf(firstline, "%d %d %d %lf", &n_graph_, &m_graph_, &nbr_demand_graph_, &scale_factor);
  if(p < 3)     {
    printf("\nError: the input file is not in the required format (1).\n");
    exit(1);
  }
  else if(p < 4) scale_factor = 1.0; 
    
  /* Allocate memory */
  or_graph_ = VECT(HighsInt, m_graph_);
  ex_graph_ = VECT(HighsInt, m_graph_);
  capac_graph_ = VECT(double, m_graph_);
  cout_graph_ = VECT(double, m_graph_);
  dem_graph_ = VECT(McnfDemand, nbr_demand_graph_);
    
  if(aggregated == 2)    {
    source_nodes_ = VECT(HighsInt, n_graph_);
    source_nodes_len_ = VECT(HighsInt, n_graph_);
    source_nodes_commodities_ = VECT(HighsInt, nbr_demand_graph_);
    source_nodes_n_ = 0;
  }
    
  /* Initialize the scale factor */
  max_cap_ = -1.0;
  max_cost_ = -1.0;

  /* Read information about arcs */
  for (i = 0; i < m_graph_; i++)     {
    ret = fgets(firstline, 1000, input_file);
    p = sscanf(firstline,  "%d %d %lf %lf", or_graph_ + i, ex_graph_ + i, cout_graph_ + i, capac_graph_ + i);
    if(p < 3)        {
      printf("\nError: the input file is not in the required format (1).\n");
      exit(1);
    }
    else if(p < 4)       {
      capac_graph_[i] = cout_graph_[i];
      cout_graph_[i] = 100.0;
    }
        
    /* Check if it is the largest cost */
    if(max_cost_ < cout_graph_[i]) max_cost_ = cout_graph_[i]; 
        
    /* Check if it is the largest capacity */
    if(capac_graph_[i] > max_cap_) max_cap_ = capac_graph_[i];

    if (zero == 0)         {
      (or_graph_[i])--;
      (ex_graph_[i])--;
    }
  }
    
  /* "Round up" max_cap */
  aux_double = max_cap_;
  max_cap_ = 1.0;
  if(aggregated != 0) while(max_cap_ < aux_double) max_cap_ *= 10.0;

  /* SCALE THE CAPACITY OF ARCS */
  if(max_cap_ > 1.0) for (i = 0; i < m_graph_; i++) capac_graph_[i] /= max_cap_; 
    
  /* "Round up" max_cost */
  aux_double = max_cost_;
  max_cost_ = 1.0;
  while(max_cost_ < aux_double) max_cost_ *= 10.0;

  /* SCALE THE LINEAR COST */
  if(max_cost_ > 1.0) for (i = 0; i < m_graph_; i++) cout_graph_[i] /= max_cost_;
    
  /* Read information about commodities */
  for (i = 0; i < nbr_demand_graph_; i++)     {
    ret = fgets(firstline, 1000, input_file);
    p = sscanf(firstline,  "%d %d %lf", &((dem_graph_ + i)->s), &((dem_graph_ + i)->t), &flot);
    if(p < 3)       {
      printf("\nError: the input file is not in the required format (2).\n");
      exit(1);
    }

    /* Apply the scaling factor read from the input file */
    (dem_graph_ + i)->val = scale_factor * flot;
        
    /* SCALE THE DEMAND OF COMMODITIES */
    if(max_cap_ > 1.0) (dem_graph_ + i)->val /= max_cap_; 
        
    if (zero == 0)         {
      ((dem_graph_ + i)->s)--;
      ((dem_graph_ + i)->t)--;
    }

    if(aggregated == 2)      {
      /* Update the number of commodities associated to the corresponding source node */
      if(source_nodes_len_[(dem_graph_ + i)->s] == 0) source_nodes_n_++;
      source_nodes_len_[(dem_graph_ + i)->s]++;
    }
  }
    
  if(aggregated == 2)    {
    /* Assign comodities to their corresponding source nodes */
    if(source_nodes_n_ > 0)      {
      /* Set the header of each group */
      source_nodes_[0] = 0;
      for (i = 1; i < n_graph_; i++)             {
	source_nodes_[i] = source_nodes_[i-1] + source_nodes_len_[i-1];
	source_nodes_len_[i-1] = 0;
      }
      source_nodes_len_[n_graph_-1] = 0;
            
      /* Set commodities in their groups */
      for (i = 0; i < nbr_demand_graph_; i++)             {
	source_nodes_commodities_[source_nodes_[(dem_graph_ + i)->s] + source_nodes_len_[(dem_graph_ + i)->s] ] = i;
	source_nodes_len_[(dem_graph_ + i)->s]++;
      }
    }
  }
    
  /* Close the input file */    
  fclose(input_file);
}

/*---------------------------- adjacency_list ----------------------------*/
void McnfOracle::adjacency_list(McnfGraph * N) {
  /* calcule la liste adjacence la (arcs adjacents a un sommet)  du graph */
  HighsInt    i, j;
  HighsInt    nb = 0;

  N->lp[0] = 0;

  if (N->oriented == 1) {
    for (i = 0; i < N->n; i++) {
      for (j = 0; j < N->m; j++) {
        if (N->org[j] == i) {
          N->la[nb] = j;
          nb++;
        }
      }
      N->lp[i + 1] = nb;
    }
  }
  else {
    for (i = 0; i < N->n; i++) {
      for (j = 0; j < N->m; j++) {
        if (N->org[j] == i || N->ex[j] == i) {
          N->la[nb] = j;
          nb++;
        }
      }
      N->lp[i + 1] = nb;
    }
  }
}

/*---------------------------- succesor_list ----------------------------*/
void McnfOracle::succesor_list(McnfGraph * N) {
  /* calcule la liste adjacence ls (sommets adjacents a un sommet)  du graph */
  HighsInt    i, k, u;

  for (k = 0; k < N->n; k++) {
    for (i = N->lp[k]; i < N->lp[k + 1]; i++) {
      u = N->la[i];
      if ((N->ls[i] = N->org[u]) == k)
        N->ls[i] = N->ex[u];
    }
  }
}

void McnfOracle::set_network_Dijkstra(McnfGraph *G, McnfDemand **dem, HighsInt *nb_dem, HighsInt oriented) {
  G->oriented = oriented;
  G->org = or_graph_;
  G->ex = ex_graph_;
  G->n = n_graph_;
  G->m = m_graph_;
  G->cost = cout_graph_;
  G->K= capac_graph_;
  G->la = VECT(HighsInt, 2 * m_graph_);
  G->ls = VECT(HighsInt, 2 * m_graph_);
  G->lp = VECT(HighsInt, n_graph_ + 1);
  adjacency_list(G);
  succesor_list(G); 

  *dem = dem_graph_;
  *nb_dem = nbr_demand_graph_;
}

/* Initialize the data structure used by the Dijkstra algorithm */
void McnfOracle::init_Dijkstra (McnfGraph * G, McnfDemand * dem, HighsInt nb_dem) {
  network_ = G;
  demand_ = dem;
  nb_demand_ = nb_dem;
  net_cost_ = VECT (double, G->m);
  prec_tableau_ = VECT (HighsInt, G->n * G->n);
  d_tableau_ = VECT (double, G->n * G->n);
  dijkstra_from_node_ = VECT (HighsInt, G->n);
  active_set_ = NULL;
  if(USE_ACTIVE_SET_) active_set_ = VECT (HighsInt, G->m);
  active_set_n_ = 0; 
}

/*---------------------------- Dijkstra ----------------------------*/
/* Find the shortest path in a graph using Dijkstra algorithm.      */
/* The adjacency lists H->la and H-​>ls are used.                    */
/*------------------------------------------------------------------*/
void ErrChem (HighsInt is, HighsInt it) {
  printf ("There is no path between %i and %i\n", is, it);
  exit (1);
}

void McnfOracle::Reset_Dijkstra (McnfGraph * H) {
  HighsInt i;
  for (i = 0; i < H->n; i++)
    dijkstra_from_node_[i] = 0;
}

void McnfOracle::Dijkstra_solve (McnfGraph * H, Vertex org, Vertex ex, double *cost) {
  HighsInt i, j, k;
  Vertex is, it;
  Arc u;
  double z;

  HighsInt *visited;
    
  is = org;
  it = ex;
  HighsInt *prec = prec_tableau_ + is * H->n;
  double *d = d_tableau_ + is * H->n;
  if (!dijkstra_from_node_[is])  {
    dijkstra_from_node_[is] = 1;

    visited = VECT (HighsInt, H->n);

    for (i = 0; i < H->n; i++)    {
      d[i] = DBL_MAX;
      prec[i] = INT_MIN;
      visited[i] = 0;
    }

    d[is] = 0;
    for (k = H->lp[is]; k < H->lp[is + 1]; k++)    {
      u = H->la[k];
      j = H->ls[k];
      prec[j] = u;
      d[j] = cost[u];
    }
    visited[is] = 1;
    i = is;           /* d(is)=0, S={is}, i=is */

    do {z = DBL_MAX; i = INT_MIN;
      for (j = 0; j < H->n; j++) {
	if (visited[j] == 0 && d[j] <= z) {
          i = j;
          z = d[j];
        }
      }
      if (i < 0) break;

      for (k = H->lp[i]; k < H->lp[i + 1]; k++)        {
	u = H->la[k];
	j = H->ls[k];
	if (visited[j] == 0)        {
          z = d[i] + cost[u];
          if (z < d[j])            {
	    d[j] = z;
	    prec[j] = u;
	  }
        }
      }
      visited[i] = 1;
    }
    while (1);
      
    if (!visited[it]) ErrChem (is, it);
    
    free(visited);
  }
  //std::cout << "finished solving Dijkstra\n";
}

/* Generate a *sparse* column using the flow of a shortest path found by the Dijkstra algorithm */
/* It is used in the disaggregated oracle */
HighsInt McnfOracle::get_solution_from_Dijkstra(McnfGraph * G, HighsPdcgm_SMatrix_CW * M, HighsInt orig, HighsInt extr, double val) {
  HighsInt u, som, *prec;

  som = extr;
  prec = prec_tableau_ + orig * G->n;
  if (prec[extr] == INT_MIN)     {
    /* No path was found by the Dijkstra algorithm */
    printf("\n### Error: No path was found by the Dijkstra algorithm.\n");
    fflush(stdout);
    return 0;
  }
    
  while (som != orig)    {
    u = prec[som];
    M->m_rwnmbs.push_back( u );
    //printf("-val=%lf\n",-val);
    M->m_coeff.push_back( -val);
    //printf("-val=%lf\n",-val);
    (M->m_nz)++;
    if (som == G->ex[u]) som = G->org[u];
    else som = G->ex[u];
  }
  //printf("In Dijkstra: \n");
  //(*M).print_SMatrix_CW();
    
  return 1;
}

/* Add to a dense column the flow of a shortest path obtained by the Dijkstra algorithm */
/* It is used in the aggregated oracle */
HighsInt McnfOracle::add_solution_from_Dijkstra(McnfGraph * G, double *vector, HighsInt orig, HighsInt extr, double val) {
  HighsInt u, som, *prec;

  som = extr;
  prec = prec_tableau_ + orig * G->n;
  if (prec[extr] == INT_MIN)    {
    /* No path was found by the Dijkstra algorithm */
    printf("\n### Error: No path was found by the Dijkstra algorithm.\n");
    fflush(stdout);
    return 0;
  }
    
  while (som != orig)    {
    u = prec[som];
    vector[u] -= val;
        
    if (som == G->ex[u]) som = G->org[u];
    else som = G->ex[u];
  }
    
  return 1;
}

