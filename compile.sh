highs_root="$HOME/HiGHS"

#highs_build_dir="${highs_root}/test"; cpp_opt_deb_flag="-O3"

highs_build_dir="${highs_root}/build"; cpp_opt_deb_flag="-g3"

source="pdcgm.cpp McnfOracle.cpp"
#pdcgm_SMatrix.cpp pdcgm_env.cpp

g++ ${cpp_opt_deb_flag} -std=c++11 -I "${highs_root}/app/"  -I "${highs_build_dir}" -I "${highs_root}/src/" -L "${highs_build_dir}/lib/" ${source} -lhighs

