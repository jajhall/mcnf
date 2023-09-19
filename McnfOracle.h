#ifndef MCNFORACLE_H
#define MCNFORACLE_H

#include "lp_data/HighsPdcgm.h"

struct McnfGraph {
    HighsInt    n;
    HighsInt    m;
    HighsInt    oriented;
    HighsInt   *org;
    HighsInt   *ex;
    HighsInt   *la;
    HighsInt   *lp;
    HighsInt   *ls;
    double *K;
    double *cost;
};

struct McnfDemand {
    HighsInt s; /* source of the demand */
    HighsInt t; /* destination of the demand */
    double val; /* quantite of flow */
};

typedef HighsInt Vertex;
typedef HighsInt Arc;
typedef double Flow;

#define EdgeEq(N,i,j)     ((N.org[i] == N.org[j] && N.ex[i] == N.ex[j])  ||   \
                 (N.org[i] == N.ex[j] && N.ex[i] == N.org[j]) )

#define EdgeEqp(N,i,j)    EdgeEq((*N),i,j)

#define max_edge_graph(n) ((n)*((n)-1)/2)

#define VECT(type, dim)  (type*) calloc((unsigned) dim, (unsigned) sizeof(type))

#define min(a,b)       ((a<b) ? a : b)
#define max(a,b)       ((a>b) ? a : b)

#define show(expr)     printf(#expr " = %g\n", expr)

class McnfOracle {
 public:

  McnfGraph *network_;
  McnfDemand *demand_;
  HighsInt nb_demand_;
  double *d_tableau_;
  double *net_cost_;
  HighsInt *prec_tableau_;
  HighsInt *dijkstra_from_node_;

  HighsInt m_graph_; 
  HighsInt n_graph_; 
  HighsInt *or_graph_; 
  HighsInt *ex_graph_;
  HighsInt nbr_demand_graph_;

  double *capac_graph_;
  double *cout_graph_;

  McnfDemand *dem_graph_;

  void adjacency_list(McnfGraph* fred);
  void succesor_list(McnfGraph* fred);

  HighsInt *active_set_;
  HighsInt active_set_n_;
  HighsInt USE_ACTIVE_SET_;

  HighsInt *source_nodes_;
  HighsInt *source_nodes_len_;
  HighsInt *source_nodes_commodities_;
  HighsInt source_nodes_n_;

  double max_cap_;
  double max_cost_;

  /* Parse the input given in the text line */
  void parseInput(int argc, char *argv[], HighsInt *aggregated, HighsInt *oriented, 
		  char DefltFlnme[], HighsInt *begin_zero, HighsInt *with_cost);

  /* Read the instance data from an input file */
  void read_instance(char file_name[], HighsInt zero, HighsInt aggregated);

  /* Header functions of the Dijkstra algorithm */
  void set_network_Dijkstra(McnfGraph *G, McnfDemand **dem, HighsInt *nb_dem, HighsInt oriented);
  void init_Dijkstra (McnfGraph * G, McnfDemand * dem, HighsInt nb_dem);
  void Reset_Dijkstra (McnfGraph * H);
  void Dijkstra_solve (McnfGraph * H, Vertex org, Vertex ex, double *cost);

  /* Generate a *sparse* column using the flow of a shortest path found by the Dijkstra algorithm */
  /* It is used in the disaggregated oracle */
  HighsInt get_solution_from_Dijkstra(McnfGraph * G, HighsPdcgm_SMatrix_CW * M, HighsInt orig, HighsInt extr, double val);

  /* Add to a dense column the flow of a shortest path obtained by the Dijkstra algorithm */
  /* It is used in the aggregated oracle */
  HighsInt add_solution_from_Dijkstra(McnfGraph * G, double *vector, HighsInt orig, HighsInt extr, double val);
};

#endif
