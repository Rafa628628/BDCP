#include <string>
using namespace std;

#ifndef _CCP_H__
#define _CCP_H__

int sign(int a);
void free_matrix(double **matrix, int r, int c);
void free_int_matrix(int **matrix, int r, int c);
void free_3d_int_matrix(int ***arr3D, int r, int c);
void free_3d_matrix(double ***arr3D, int r, int c);
double **alloc_matrix(int r, int c);
double ***alloc_3d_matrix(int r, int c, int h);
int **alloc_int_matrix(int r, int c);
int ***alloc_3d_int_matrix(int r, int c, int h);
double Euclidean_distance(double *p, int d);
double det(double *a, double *b);
void line_intersection(double *intersection, int *need_flag, double **line1, double **line2);
//double ** cal_distance(double **time_series, char norm, int len, int d);
double get_auto_corr(double **time_series, int time_begin, int time_end, int d, int order);
class cp
{
	public:
	    double **time_series;
		int time_series_len; 
		int time_series_d;
		int min_segment_size;
		double pt;
		bool is_Euclidean_distance;
		double **dist;
        int block_size;
        int * change_point_set;
        int * change_point_num;
		
		cp(double *my_time_series, int my_time_series_len, int my_time_series_d, int my_min_segment_size, double my_pt, bool my_is_Euclidean_distance, double *my_dist, int my_block_size, int * my_change_point_set, int * my_change_point_num);
		double calculate_V(int * series_index, int time_begin, int m, int l);
		void cp_origin();
		double cal_threshold_cp(int time_begin, int time_end);
		void golden_section_search(int *best_m_l, double *best_v, int *p, int *u, int time_begin, int time_end, int * series_index, int tol);
		void cp_gss();
		double cal_threshold_cp_gss(int time_begin, int time_end);
		void cp_powell();
		double cal_threshold_cp_powell(int time_begin, int time_end);
		int get_block_size(int time_begin, int time_end);
		int get_qt(double **my_time_series, int time_begin, int time_end);
		void powell(int *best_m_l, double *best_v, int time_begin, int time_end, int * series_index);
		void cal_intersections(int *a, int *b, int *p, int *u, int time_begin, int time_end);		
};
void run(double *my_time_series, int my_time_series_len, int my_time_series_d, int my_min_segment_size, double my_pt, bool my_is_Euclidean_distance, double *my_dist, int my_block_size, int * my_change_point_set, int * my_change_point_num, char method);

//
//void init(double **my_time_series, int my_time_serie_len, int my_time_series_d, int my_min_segment_size, double my_pt, char my_norm);
//void cp_origin(int *change_point_set, int *change_point_num);
//void cp_gss(int *change_point_set, int *change_point_num);
//void cp_powell(int *change_point_set, int *change_point_num);

#endif	// CCP_H__ 
