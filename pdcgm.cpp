#include "Highs.h"
#include "McnfOracle.h"

const bool dev_run = true;

static void userColOracle(const int col_oracle_type, const char* message,
                             const HighsColOracleDataOut* data_out,
                             HighsColOracleDataIn* data_in,
                             void* user_col_oracle_data);
int main(int argc, char *argv[]) {
  // Auxiliary variables
  HighsInt nb_dem, aggregated, oriented, begin_zero, with_cost;

  // Graph representing the network
  McnfGraph G; 
  // Demand of each commodity
  McnfDemand *dem; 

  // Instance file name
  char DefltFlnme[1000];
  // Output file pointer
  FILE *out_file;

  // Parse the command line receive as an input
  parseInput(argc, argv, &aggregated, &oriented, DefltFlnme, &begin_zero, &with_cost);

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
