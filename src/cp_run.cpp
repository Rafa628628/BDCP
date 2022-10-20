# include <Rcpp.h>
# include "ccp.h"

# include<iostream>
# include<stdlib.h>
# include<cmath>
# include<ctime>
# include <string>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector cp_run(NumericVector my_time_series_raw,
                     int my_time_series_len,
                     int my_time_series_d,
                     int my_min_segment_size,
                     double my_pt,
                     bool my_is_Euclidean_distance,
                     NumericVector my_dist_raw,
                     int my_block_size,
                     NumericVector my_change_point_set_raw,
                     NumericVector my_change_point_num_raw,
                     char method)

{
  std::vector<double> my_time_series_temp;
  std::copy(my_time_series_raw.begin(), my_time_series_raw.end(), std::back_inserter(my_time_series_temp));
  double* my_time_series = &my_time_series_temp[0];

  std::vector<double> my_dist_temp;
  std::copy(my_dist_raw.begin(), my_dist_raw.end(), std::back_inserter(my_dist_temp));
  double* my_dist = &my_dist_temp[0];

  std::vector<int> my_change_point_set_temp;
  std::copy(my_change_point_set_raw.begin(), my_change_point_set_raw.end(), std::back_inserter(my_change_point_set_temp));
  int* my_change_point_set = &my_change_point_set_temp[0];

  std::vector<int> my_change_point_num_temp;
  std::copy(my_change_point_num_raw.begin(), my_change_point_num_raw.end(), std::back_inserter(my_change_point_num_temp));
  int* my_change_point_num = &my_change_point_num_temp[0];

  run(my_time_series, my_time_series_len, my_time_series_d, my_min_segment_size, my_pt, my_is_Euclidean_distance, my_dist, my_block_size, my_change_point_set, my_change_point_num, method);

  NumericVector res = NumericVector(my_change_point_set, my_change_point_set + my_change_point_num[0]+1);

  my_dist_temp.clear();
  my_time_series_temp.clear();
  my_change_point_set_temp.clear();
  my_change_point_num_temp.clear();

  return res;
}

