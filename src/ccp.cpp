#include "ccp.h"
#include "stdlib.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
# define PI 3.1415926
using namespace std;

int sign(int a)
{
	if(a>0)
	{
		return 1;
	}
	else if(a<0)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

void free_matrix(double **matrix, int r, int c) {
    /* free a matrix with r rows and c columns */
    int i;
    for (i = 0; i < r; i++) {
//        cout<<i<<endl;
        free(matrix[i]);
    }
    free(matrix);
}

void free_int_matrix(int **matrix, int r, int c) {
    /* free a matrix with r rows and c columns */
    int i;
    for (i = 0; i < r; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void free_3d_int_matrix(int ***arr3D, int r, int c) {
    int i, j;

    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            free(arr3D[i][j]);
        }
        free(arr3D[i]);
    }
    free(arr3D);
}

void free_3d_matrix(double ***arr3D, int r, int c) {
    int i, j;

    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            free(arr3D[i][j]);
        }
        free(arr3D[i]);
    }
    free(arr3D);
}

double **alloc_matrix(int r, int c) {
    /* allocate a matrix with r rows and c columns */
    int i;
    double **matrix;
    matrix = (double **) calloc((size_t) r, sizeof(double *));
    for (i = 0; i < r; i++)
        matrix[i] = (double *) calloc((size_t) c, sizeof(double));
    return matrix;
}

double ***alloc_3d_matrix(int r, int c, int h) {
    /* allocate a 3D matrix with r rows, c columns, h levels */
    double ***arr3D;
    int i, j;

    arr3D = (double ***) malloc(r * sizeof(double **));

    for (i = 0; i < r; i++) {
        arr3D[i] = (double **) malloc(c * sizeof(double *));
        for (j = 0; j < c; j++) {
            arr3D[i][j] = (double *) malloc(h * sizeof(double));
        }
    }
    return arr3D;
}

int **alloc_int_matrix(int r, int c) {
    /* allocate a matrix with r rows and c columns */
    int i;
    int **matrix;
    matrix = (int **) calloc((size_t) r, sizeof(int *));
    for (i = 0; i < r; i++)
        matrix[i] = (int *) calloc((size_t) c, sizeof(int));
    return matrix;
}

int ***alloc_3d_int_matrix(int r, int c, int h) {
    /* allocate a 3D matrix with r rows, c columns, h levels */
    int ***arr3D;
    int i, j;

    arr3D = (int ***) malloc(r * sizeof(int **));

    for (i = 0; i < r; i++) {
        arr3D[i] = (int **) malloc(c * sizeof(int *));
        for (j = 0; j < c; j++) {
            arr3D[i][j] = (int *) malloc(h * sizeof(int));
        }
    }
    return arr3D;
}


double Euclidean_distance(double *p, int d)
{
	double temp=0;
	int i;
	for(i=0;i<d;i++)
	{
		temp += pow(p[i], 2);
	}
	temp=pow(temp/d, 0.5);
	return temp;
}


double det(double *a, double *b)
{
	return a[0]*b[1] - a[1] * b[0];
}

void line_intersection(double *intersection, int *need_flag, double **line1, double **line2)
{
	double *xdiff, *ydiff, *d;
	double div;
	xdiff=(double*) malloc(2*sizeof(double));
	ydiff=(double*) malloc(2*sizeof(double));
	d=(double*) malloc(2*sizeof(double));
	
	xdiff[0] = line1[0][0] - line1[1][0];
	xdiff[1] = line2[0][0] - line2[1][0];
	ydiff[0] = line1[0][1] - line1[1][1];
	ydiff[1] = line2[0][1] - line2[1][1];

    div = det(xdiff, ydiff);
    if(div == 0)
    {
    	need_flag[0] = 0;
    	free(xdiff);
    	free(ydiff);
    	free(d);
		return;  	
	}
	else
	{
		d[0] = det(line1[0], line1[1]);
		d[1] = det(line2[0], line2[1]);

	    intersection[0] = det(d, xdiff) / div;
	    intersection[1] = det(d, ydiff) / div;
	    need_flag[0] = 1;
	    free(xdiff);
    	free(ydiff);
    	free(d);
	    return; 
	}
}

//double ** cal_distance(double **time_series, char norm, int len, int d)
//{
//	double *temp;
//	temp = (double*) malloc(d*sizeof(double));
//	int i,j,k;
//	double ** dist;
//	dist = alloc_matrix(len, len);
//	for (i=0;i<len;i++)
//	{
//		for (j=i;j<len;j++)
//		{
//			for(k=0;k<d;k++)
//			{
//				temp[k] = time_series[i][k] - time_series[j][k];
//			}
//			if(norm == '1')
//			{
//				dist[i][j] = norm_1(temp, d);
//				dist[j][i] = norm_1(temp, d);
//			}
//			else
//			{
//				dist[i][j] = norm_2(temp, d);
//				dist[j][i] = norm_2(temp, d);
//			}
//		}
//	}
//	free(temp);
//	return dist;
//}

double get_auto_corr(double **time_series, int time_begin, int time_end, int d, int order)
{
	int i;
	double *time_series_mean;
	time_series_mean = (double*) malloc(d*sizeof(double));
	double time_series_var=0;
	double temp_mean=0;
	double auto_corr=0;
	double temp_auto_corr = 0.0001;
	
	double *temp;
	temp = (double*) malloc(d*sizeof(double));
	int j;
	
	for(j=0;j<d;j++)
	{
		for(i=time_begin;i<time_end;i++)
		{
			temp_mean += time_series[i][j];
		}
		time_series_mean[j] = temp_mean / (time_end - time_begin);
	}
	for(i=time_begin;i<time_end;i++)
	{	
		for(j=0;j<d;j++)
		{
			temp[j] = time_series[i][j] - time_series_mean[j];
		}
        time_series_var += pow(Euclidean_distance(temp, d), 2);
	}
	for(i=time_begin;i<time_end-order;i++)
	{	
		temp_auto_corr = 0;
		for(j=0;j<d;j++)
		{
			temp_auto_corr += (time_series[i][j] - time_series_mean[j])*(time_series[i+order][j] - time_series_mean[j]);
		}		
		auto_corr += temp_auto_corr / time_series_var;
	}
	free(time_series_mean);
	free(temp);
	return auto_corr;	
}

//class cp
//{
//	public:
//		int time_series_len;
//		int time_series_d;
//		double **time_series;
//		int min_segment_size;
//		double pt;
//		double **dist;
//		char norm;
//
//		void init(double **my_time_series, int my_time_serie_len, int my_time_series_d, int my_min_segment_size, double my_pt, char my_norm);
//		double calculate_V(double **mydist, int time_begin, int m, int l);
//		void cp_origin(int *change_point_set, int *change_point_num);
//		double cal_threshold_cp(double **time_series, int time_begin, int time_end);
//		void golden_section_search(int *best_m_l, double *best_v, int *p, int *u, int time_begin, int time_end, double **dist, int tol);
//		void cp_gss(int *change_point_set, int *change_point_num);
//		double cal_threshold_cp_gss(double **time_series, int time_begin, int time_end);
//		void cp_powell(int *change_point_set, int *change_point_num);
//		double cal_threshold_cp_powell(double **time_series, int time_begin, int time_end);
//		int get_block_size(double **time_series, int time_begin, int time_end);
//		int get_qt(double **time_series, int time_begin, int time_end);
//		void powell(int *best_m_l, double *best_v, int time_begin, int time_end, double **dist);
//		void cal_intersections(int *a, int *b, int *p, int *u, int time_begin, int time_end);
//};

cp::cp(double *my_time_series, int my_time_series_len, int my_time_series_d, int my_min_segment_size, double my_pt, bool my_is_Euclidean_distance, double *my_dist, int my_block_size, int * my_change_point_set, int * my_change_point_num)
{

//    cout<<"init"<<endl;
	int i,j;
	double **ts_copy;
	ts_copy = alloc_matrix(my_time_series_len, my_time_series_d);
	for(i=0;i<my_time_series_len;i++)
	{
		for(j=0;j<my_time_series_d;j++)
		{
			ts_copy[i][j]=my_time_series[i*my_time_series_d + j];
		}
	}
	time_series = ts_copy;
//    cout<<"time_series"<<endl;
//    for(i=0;i<my_time_serie_len;i++)
//	{
//		for(j=0;j<my_time_series_d;j++)
//		{
//			cout<<time_series[i][j]<<' ';
//		}
//		cout<<endl;
//	}
//	cout<<endl;
	time_series_len = my_time_series_len;
	time_series_d = my_time_series_d;
	min_segment_size = my_min_segment_size;
	pt = my_pt;
	is_Euclidean_distance = my_is_Euclidean_distance;

	double **dist_copy;
	dist_copy = alloc_matrix(my_time_series_len, my_time_series_len);
	for(i=0;i<my_time_series_len;i++)
	{
		for(j=0;j<my_time_series_len;j++)
		{
			dist_copy[i][j]=my_dist[i * my_time_series_len + j];
		}
	}
	dist = dist_copy;
	block_size = my_block_size;
	change_point_set = my_change_point_set;
	change_point_num = my_change_point_num;
}

double cp::calculate_V(int *series_index, int time_begin, int m, int l)
{
	int x,y,i,j;
	double temp,c_1,c_2,v;
	
	if(m == l)
	{
		return 0;
	}
	temp=0;
	for(i=time_begin;i<l;i++)
	{
		for(j=time_begin;j<l;j++)
		{
			c_1=0;
			c_2=0;
			for(x=time_begin;x<m;x++)
			{	
				if(dist[series_index[i]][series_index[x]]<=dist[series_index[i]][series_index[j]])
				{
					c_1++;
				}
			}
				
			for(y=m;y<l;y++)
			{
				if(dist[series_index[i]][series_index[y]]<=dist[series_index[i]][series_index[j]])
				{
					c_2++;
				}
			}
			temp = temp + pow((c_1/(m-time_begin) - c_2/(l-m)),2);
			
		}
	}
	v = temp * (m - time_begin) * (l - m) / pow((l - time_begin), 3);
	return v;
		
}

void cp::cal_intersections(int *a, int *b, int *p, int *u, int time_begin, int time_end)
{
	double **line0, ***line_set, **intersections;
	int **need_flag;
	int i,j;
	
	line0= alloc_matrix(2, 2);
	line_set=alloc_3d_matrix(4, 2, 2);
	need_flag=alloc_int_matrix(4,1);
	intersections=alloc_matrix(4, 2);
	
	line0[0][0] = double(p[0]);
	line0[0][1] = double(p[1]);
	line0[1][0] = double(p[0] + u[0]);
	line0[1][1] = double(p[1] + u[1]);
	
	line_set[0][0][0] = double(0);
	line_set[0][0][1] = double(time_end);
	line_set[0][1][0] = double(time_end); 
	line_set[0][1][1] = double(time_end);
	
	line_set[1][0][0] = double(time_begin + min_segment_size);
	line_set[1][0][1] = double(0);
	line_set[1][1][0] = double(time_begin + min_segment_size); 
	line_set[1][1][1] = double(time_end);
	
	line_set[2][0][0] = double(time_end - min_segment_size);
	line_set[2][0][1] = double(0);
	line_set[2][1][0] = double(time_end - min_segment_size); 
	line_set[2][1][1] = double(time_end);
	
	line_set[3][0][0] = double(0);
	line_set[3][0][1] = double(0);
	line_set[3][1][0] = double(time_end); 
	line_set[3][1][1] = double(time_end);
	
	for(i=0;i<4;i++)
	{
		line_intersection(intersections[i], need_flag[i], line0, line_set[i]);
	}
	
	for(i=0;i<4;i++)
	{
		if(need_flag[i][0] == 1)
		{
			if((intersections[i][0] < time_begin+min_segment_size) | (intersections[i][0] > time_end-min_segment_size) | (intersections[i][1] < intersections[i][0]) | (intersections[i][1] > time_end))
			{
				need_flag[i][0] = 0;
			}
		}
	}
	
	for(i=0;i<4;i++)
	{
		if(need_flag[i][0] == 1)
		{
			for(j=i+1;j<4;j++)
			{
				if(need_flag[j][0] == 1)
				{
					if(intersections[i][0] == intersections[j][0] && intersections[i][1] == intersections[j][1])
					{
						need_flag[j][0] = 0;
					}
				}
			}
			
		}
		
	}
	
	j=0;
	for(i=0;i<4;i++)
	{
		if(need_flag[i][0] == 1)
		{
		    if(j==3)
		    {
		        j +=1;
		    }
            else if(j==2)
            {
                j +=1;
            }
			else if(j==1)
			{
				b[0] = int(intersections[i][0]);
				b[1] = int(intersections[i][1]);
                j+=1;
			}
			else if(j==0)
			{
				a[0] = int(intersections[i][0]);
				a[1] = int(intersections[i][1]);
				j +=1;
			}
		}
	}
	
	if(j != 2)
	{
		cout<<"---------------------------"<<endl;
		cout<<"inetrsection numbers wrong"<<j<<endl;
		cout<<"p"<<p[0]<<","<<p[1]<<endl;
		cout<<"u"<<u[0]<<","<<u[1]<<endl;
		cout<<"intersections[0]"<<intersections[0][0]<<","<<intersections[0][1]<<endl;
		cout<<"intersections[1]"<<intersections[1][0]<<","<<intersections[1][1]<<endl;
		cout<<"intersections[2]"<<intersections[2][0]<<","<<intersections[2][1]<<endl;
		cout<<"intersections[3]"<<intersections[3][0]<<","<<intersections[3][1]<<endl;
		cout<<"need_flag[0]"<<need_flag[0][0]<<endl;
		cout<<"need_flag[1]"<<need_flag[1][0]<<endl;
		cout<<"need_flag[2]"<<need_flag[2][0]<<endl;
		cout<<"need_flag[3]"<<need_flag[3][0]<<endl;
		cout<<"time_begin: "<<time_begin<<endl;
		cout<<"time_end: "<<time_end<<endl;
	}
//    cout<<"intersection end"<<endl;
	free_matrix(line0,2,2);
	free_matrix(intersections,4,2);
	free_int_matrix(need_flag,4,1);
	free_3d_matrix(line_set,4,2);
//	cout<<"intersection end 2"<<endl;
	return;
}

void cp::golden_section_search(int *best_m_l, double *best_v, int *p, int *u, int time_begin, int time_end, int *series_index, int tol)
{
//    cout<<"gss"<<endl;
    if((u[0]==0) && (u[1]==0))
    {
        best_m_l[0] = p[0];
        best_m_l[1] = p[1];
        best_v[0] = calculate_V(series_index, time_begin, p[0], p[1]);
        return;
    }
	double invphi, invphi2, yc, yd;
	int h_flag, *a, *b, *h, *c, *d, i;
	
	invphi = (pow(5,0.5) - 1.0) / 2.0;
	invphi2 = (3.0 - pow(5,0.5)) / 2.0;
	a=(int*) malloc(2*sizeof(int));
	b=(int*) malloc(2*sizeof(int));
	c=(int*) malloc(2*sizeof(int));
	d=(int*) malloc(2*sizeof(int));
	h=(int*) malloc(2*sizeof(int));
//	cout<<"cal_intersection"<<endl;
	cal_intersections(a,b,p,u,time_begin,time_end);
//	cout<<"cal_intersection_end"<<endl;
	
	h[0] = b[0] - a[0];
	h[1] = b[1] - a[1];
	
	h_flag = max(abs(h[0]), abs(h[1]));
	
	c[0] = a[0] + sign(h[0]) * int(abs(invphi2 * h[0]));
	c[1] = a[1] + sign(h[1]) * int(abs(invphi2 * h[1]));
	d[0] = a[0] + sign(h[0]) * ceil(abs(invphi * h[0]));
	d[1] = a[1] + sign(h[1]) * ceil(abs(invphi * h[1]));
	
	yc = calculate_V(series_index, time_begin, c[0], c[1]);
	yd = calculate_V(series_index, time_begin, d[0], d[1]);
	
//	cout<<"p: "<<p[0]<<","<<p[1]<<endl;
//	cout<<"u: "<<u[0]<<","<<u[1]<<endl;
//	cout<<"a: "<<a[0]<<","<<a[1]<<endl;
//	cout<<"b: "<<b[0]<<","<<b[1]<<endl;
//	cout<<"c: "<<c[0]<<","<<c[1]<<endl;
//	cout<<"d: "<<d[0]<<","<<d[1]<<endl;
	
	if(h_flag <= 1)
	{
		if(yc > yd)
		{
			best_m_l[0] = c[0];
			best_m_l[1] = c[1];
			best_v[0] = yc;
            free(a);
            free(b);
            free(c);
            free(d);
            free(h);
			return;
		}
		else
		{
			best_m_l[0] = d[0];
			best_m_l[1] = d[1];
			best_v[0] = yd;
			free(a);
            free(b);
            free(c);
            free(d);
            free(h);
			return;
		}
	}
	
	while(1)
	{
//	    cout<<"while"<<endl;
        if(yc > yd)
        {
        	b[0] = d[0];
        	b[1] = d[1];
        	d[0] = c[0];
        	d[1] = c[1];
	        yd = yc;
	        h[0] = b[0] - a[0];
			h[1] = b[1] - a[1];
	        h_flag = max(abs(h[0]), abs(h[1]));
	        
	        c[0] = a[0] + sign(h[0]) * int(abs(invphi2 * h[0]));
			c[1] = a[1] + sign(h[1]) * int(abs(invphi2 * h[1]));
	        yc = calculate_V(series_index, time_begin, c[0], c[1]);
		}
    	else
		{
			a[0] = c[0];
        	a[1] = c[1];
        	c[0] = d[0];
        	c[1] = d[1];
	        yc = yd;
 			h[0] = b[0] - a[0];
			h[1] = b[1] - a[1];
	        h_flag = max(abs(h[0]), abs(h[1]));
	        
	        d[0] = a[0] + sign(h[0]) * ceil(abs(invphi * h[0]));
			d[1] = a[1] + sign(h[1]) * ceil(abs(invphi * h[1]));
			yd = calculate_V(series_index, time_begin, d[0], d[1]);
		}
        

	    if(h_flag <= tol)
		{
//		    cout<<"tol"<<endl;
		    int best_m_in_tol=0, best_l_in_tol=0, temp_m, temp_l;
		    double best_v_in_tol=-1, temp_v;
		    for(i=0;i<=abs(c[0] - d[0]);i++)
		    {
//		        cout<<abs(c[0] - d[0])<<endl;
//		        cout<<i<<endl;
		        temp_m = c[0] + sign(d[0] - c[0]) * i;
		        if(c[0] - d[0] == 0)
		        {
		            temp_l = c[1];
		        }
		        else
		        {
		            temp_l = c[1] + int(i * (d[1] - c[1]) / abs(c[0] - d[0]));
		        }
		        temp_v = calculate_V(series_index, time_begin, temp_m, temp_l);
//		        cout<<"v end"<<endl;
		        if(temp_v > best_v_in_tol)
		        {
		            best_m_in_tol = temp_m;
		            best_l_in_tol = temp_l;
		            best_v_in_tol = temp_v;
		        }
		    }
		    best_m_l[0] = best_m_in_tol;
            best_m_l[1] = best_l_in_tol;
            best_v[0] = best_v_in_tol;
            free(a);
            free(b);
            free(c);
            free(d);
            free(h);
            return;
//			if(yc > yd)
//			{
//				best_m_l[0] = c[0];
//				best_m_l[1] = c[1];
//				best_v[0] = yc;
//                free(a);
//                free(b);
//                free(c);
//                free(d);
//                free(h);
//				return;
//			}
//			else
//			{
//				best_m_l[0] = d[0];
//				best_m_l[1] = d[1];
//				best_v[0] = yd;
//				free(a);
//                free(b);
//                free(c);
//                free(d);
//                free(h);
//				return;
//			}
		}
	}	
} 

void cp::powell(int *best_m_l, double *best_v, int time_begin, int time_end, int* series_index)
{
	int **P, **U, *temp;
	int i, loop_count=0, loop_max = 20;
    double v_old, v_new;

	temp = (int*) malloc(2*sizeof(int));
	P=alloc_int_matrix(3, 2);
	U=alloc_int_matrix(2, 2);
	
	P[0][0] = time_begin + min_segment_size;
	P[0][1] = time_end;
	
	U[0][0] = 1;
	U[0][1] = 0;
	U[1][0] = 0;
	U[1][1] = 1;
	
	while(1)
	{
        loop_count +=1;
		for(i=0;i<2;i++)
		{
			golden_section_search(temp, best_v, P[i], U[i], time_begin, time_end, series_index, 2);
			v_old = calculate_V(series_index, time_begin, P[i][0], P[i][1]);
			v_new = calculate_V(series_index, time_begin, temp[0], temp[1]);
			if(v_new > v_old)
			{
			    P[i+1][0] = temp[0];
			    P[i+1][1] = temp[1];
			}
			else
			{
			    P[i+1][0] = P[i][0];
			    P[i+1][1] = P[i][1];

			}
		}
		U[0][0] = U[1][0];
		U[0][1] = U[1][1];
		U[1][0] = P[2][0] - P[0][0];
		U[1][1] = P[2][1] - P[0][1];
		
		if(U[1][0] == 0 && U[1][1] == 0)
		{
			best_m_l[0] = P[0][0];
			best_m_l[1] = P[0][1];
            free_int_matrix(P,3,2);
            free_int_matrix(U,2,2);
            free(temp);
//            cout<<"loop time: "<<loop_count<<endl;
			return;
		}
		else if(loop_count >= loop_max)
		{
        best_m_l[0] = P[2][0];
        best_m_l[1] = P[2][1];

        free_int_matrix(P,3,2);
        free_int_matrix(U,2,2);
        free(temp);
        cout<<"loop time-------------------------: "<<loop_count<<endl;
        return;
		}
		else
		{
			golden_section_search(temp, best_v, P[0], U[1], time_begin, time_end, series_index, 2);

			v_old = calculate_V(series_index, time_begin, P[2][0], P[2][1]);
			v_new = calculate_V(series_index, time_begin, temp[0], temp[1]);
			if(v_new >= v_old)
			{
			    P[0][0] = temp[0];
			    P[0][1] = temp[1];
			}
			else
			{
			    P[0][0] = P[2][0];
			    P[0][1] = P[2][1];
			}

//			cout<<"temp"<<temp[0]<<","<<temp[1]<<endl;
//			cout<<"v_new: "<<v_new<<endl;
//			cout<<"v_old: "<<v_old<<endl;
//			cout<<"2"<<endl;
//			cout<<"P[0]"<<P[0][0]<<","<<P[0][1]<<endl;
//			cout<<"P[1]"<<P[1][0]<<","<<P[1][1]<<endl;
//			cout<<"P[2]"<<P[2][0]<<","<<P[2][1]<<endl;
//			cout<<"U[0]"<<U[0][0]<<","<<U[0][1]<<endl;
//			cout<<"U[1]"<<U[1][0]<<","<<U[1][1]<<endl;
		}
	}
}

int cp::get_qt(double **my_time_series, int time_begin, int time_end)
{
	int qt_1;
	int qt_2;
	double auto_corr;
	auto_corr = get_auto_corr(my_time_series, time_begin, time_end, time_series_d, 1);
	qt_1 = (int)pow(3.0/2.0*double(time_end-time_begin), 1/3) * pow(pow((2 * auto_corr/(1-auto_corr * auto_corr)), 2), 1/3);
	qt_2 = (int)(8 * pow(double(time_end-time_begin) / 100, 1/3));
	if(qt_1>qt_2)
	{
		return qt_2;
	}		
	else
	{
		return qt_1;
	}
}

int cp::get_block_size(int time_begin, int time_end)
{
    if(!is_Euclidean_distance)
    {
        return block_size;
    }
    else
    {
        double **time_series_2;
        int i, j, qt_1, qt_2;

        time_series_2=alloc_matrix(time_series_len, time_series_d);
        for(i=time_begin;i<time_end;i++)
        {
            for(j=0;j<time_series_d;j++)
            {
                time_series_2[i][j] = pow(time_series[i][j], 2);
            }
        }

        qt_1 = get_qt(time_series, time_begin, time_end);
        qt_2 = get_qt(time_series_2, time_begin, time_end);
        if(qt_1<1 && qt_2<1)
        {
            free_matrix(time_series_2,time_series_len,time_series_d);
            return 1;
        }
        else
        {
            if(qt_1>qt_2)
            {
                free_matrix(time_series_2,time_series_len,time_series_d);
                return qt_1;
            }
            else
            {
                free_matrix(time_series_2,time_series_len,time_series_d);
                return qt_2;
            }
        }
    }



}

double cp::cal_threshold_cp(int time_begin, int time_end)
 {
//    cout<<"cal_threshold"<<endl;
	int r,b,l,i,j,m;
	int B, L, R=50;
	double *sample_best;
	int * sample_series_index;
	int sample_idx;
	double v,best;

//    cout<<"dist"<<endl;
//    for(i=0;i<time_series_len;i++)
//	{
//		for(j=0;j<time_series_len;j++)
//		{
//			cout<<dist[i][j]<<' ';
//		}
//		cout<<endl;
//	}
//	cout<<endl;

	B = get_block_size(time_begin, time_end);
	L = (time_end-time_begin)/B + 1;
//	cout<<"time_begin: "<<time_begin<<endl;
//	cout<<"time_end: "<<time_end<<endl;
//	cout<<"B, L"<<B<<','<<L<<endl;
		
	sample_best=(double*) malloc(R*sizeof(double));
	
	sample_series_index = (int*) malloc((time_end-time_begin)*sizeof(int*));

    srand(123);
	for(r=0;r<R;r++)
	{
//		cout<<"r: "<<r<<endl;
		for(l=0;l<L;l++)
		{
			sample_idx =  time_begin + rand()%(time_end-time_begin - B + 1);
//			cout<<sample_idx<<' ';
			for(b=0;b<B;b++)
			{
				if(B*l+b<time_end-time_begin)
				{
				    sample_series_index[B*l+b] = sample_idx+b;
				}
				else
				{
					break;
				}
			}
			if(B*l+b>=time_end-time_begin)
			{
				break;
			}
		}
		best = 0;
		for(m = min_segment_size; m<time_end - time_begin -min_segment_size; m++)
		{
			for(l = m + 1; l<time_end - time_begin; l++)
			{
				v = calculate_V(sample_series_index, 0, m, l);
				if(v > best)
				{
					best=v;
				}
			}
		}
		sample_best[r]=best;
	}

//	cout<<"sample_best"<<endl;
//    for(i=0;i<R;i++)
//	{
//        cout<<sample_best[i]<<' ';
//	}
//	cout<<endl;

//	cout<<"sample_time_series"<<endl;
//    for(i=0;i<time_end - time_begin;i++)
//	{
//		for(j=0;j<time_series_d;j++)
//		{
//			cout<<sample_time_series[i][j]<<' ';
//		}
//		cout<<endl;
//	}
//	cout<<endl;

	int sort_num= int(R * pt);
//	cout<<"sort_num: "<<sort_num<<endl;
	double *max_sort;
	double threshold;
	max_sort=(double*) malloc(sort_num*sizeof(double));
	for(i=0;i<sort_num;i++)
	{
	    max_sort[i] = 0;
	}

	for(r=0;r<R;r++)
	{
//	    cout<<"sample_best: "<<sample_best[r]<<endl;
		for(i=0;i<sort_num;i++)
		{
			if(sample_best[r] > max_sort[i])
			{
			    if(sort_num > 1)
			    {
			        for(j=sort_num-2;j>=i;j--)
				    {
					    max_sort[j + 1] = max_sort[j];
				    }
				    max_sort[i] = sample_best[r];
			    }
			    else
			    {
			        max_sort[0] = sample_best[r];
			    }
				break;
			}
			
		}
//		cout<<"max_sort: "<<max_sort[sort_num - 1]<<endl;
	}
//	cout<<"threshold end"<<endl;
	threshold=max_sort[sort_num - 1];
	free(max_sort);
    free(sample_best);
    free(sample_series_index);
//    cout<<"threshold end 0"<<endl;
//	free_matrix(sample_time_series, time_end-time_begin, time_series_d);
//	cout<<"threshold end 1"<<endl;
//	free_matrix(sample_dist, time_end - time_begin, time_end - time_begin);
//    cout<<"threshold end 2"<<endl;
	return threshold;
}

double cp::cal_threshold_cp_gss(int time_begin, int time_end)
{
	int r,b,l,i,j,m;
	int B, L, R=50;
	double *sample_best;
//	double **sample_time_series;
//	double **sample_dist;
    int * sample_series_index;
	int sample_idx;
	
	B = get_block_size(time_begin, time_end);
	L = (time_end-time_begin)/B + 1;
//    cout<<"B, L"<<B<<','<<L<<endl;
		
	sample_best=(double*) malloc(R*sizeof(double));
	
//	sample_time_series=alloc_matrix(time_end-time_begin, time_series_d);
    sample_series_index = (int *) malloc((time_end-time_begin) * sizeof(int));
	
	double best;
	double *sample_this_m_best_v;
	int *sample_this_m_best_m_l;

	sample_this_m_best_v =(double*) malloc(1*sizeof(double));
	sample_this_m_best_m_l =(int*) malloc(2*sizeof(int));
	
	int *p, *u;
	p =(int*) malloc(2*sizeof(int));
	u =(int*) malloc(2*sizeof(int));

    srand(123);
	for(r=0;r<R;r++)
	{
//		cout<<"r: "<<r<<endl;
//        cout<<"r";
		for(l=0;l<L;l++)
		{
			sample_idx =  time_begin + rand()%(time_end-time_begin - B + 1);
//			cout<<sample_idx<<' ';
			for(b=0;b<B;b++)
			{	
				if(B*l+b<time_end-time_begin)
				{
				    sample_series_index[B*l+b] = sample_idx + b;
				}
				else
				{
					break;
				}		
			}
			if(B*l+b>=time_end-time_begin)
			{
				break; 
			}
		}
		best = 0;
		for(m = min_segment_size; m<time_end-time_begin - min_segment_size; m++)
		{
			p[0] = m;
			p[1] = m+1;
			u[0] = 0;
			u[1] = 1;					
			golden_section_search(sample_this_m_best_m_l, sample_this_m_best_v, p, u, 0, time_end-time_begin, sample_series_index, 2);
			if(sample_this_m_best_v[0] > best)
			{
				best = sample_this_m_best_v[0];
			}		
		}
		sample_best[r]=best;
	}
//	cout<<"sample_best"<<endl;
//    for(i=0;i<R;i++)
//	{
//        cout<<sample_best[i]<<' ';
//	}
//	cout<<endl;

//	cout<<"sample_time_series"<<endl;
//    for(i=0;i<time_end - time_begin;i++)
//	{
//	    cout<<i<<":";
//		for(j=0;j<time_series_d;j++)
//		{
//			cout<<sample_time_series[i][j]<<' ';
//		}
//	}
//	cout<<endl;

	int sort_num= int(R * pt);
	double *max_sort;
	double threshold;
	max_sort=(double*) malloc(sort_num*sizeof(double));
	for(i=0;i<sort_num;i++)
	{
	    max_sort[i] = 0;
	}

	for(r=0;r<R;r++)
	{
		for(i=0;i<sort_num;i++)
		{
			if(sample_best[r] > max_sort[i])
			{
			    if(sort_num > 1)
			    {
			        for(j=sort_num-2;j>=i;j--)
				    {
					    max_sort[j + 1] = max_sort[j];
				    }
				    max_sort[i] = sample_best[r];
			    }
			    else
			    {
			        max_sort[0] = sample_best[r];
			    }
				break;
			}

		}
//		cout<<"max_sort: "<<max_sort[sort_num - 1]<<endl;
	}

//	cout<<"max_sort"<<endl;
//    for(i=0;i<sort_num;i++)
//	{
//        cout<<max_sort[i]<<' ';
//	}
//	cout<<endl;
	threshold=max_sort[sort_num - 1];
	free(p);
	free(u);
//    cout<<"threshold end 0"<<endl;
	free(sample_this_m_best_v);
	free(sample_this_m_best_m_l);
	free(max_sort);
    free(sample_best);
//    cout<<"threshold end 1"<<endl;
    free(sample_series_index);
//    cout<<"threshold end 2"<<endl;
//    cout<<"threshold end 1"<<endl;
//	free_matrix(sample_time_series, time_end-time_begin, time_series_d);
//    cout<<"threshold end 2"<<endl;
//	free_matrix(sample_dist, time_end - time_begin, time_end - time_begin);
//    cout<<"threshold end 3"<<endl;
	return threshold;
}

double cp::cal_threshold_cp_powell(int time_begin, int time_end)
{
	int r,b,l,i,j;
	int B, L, R=50;
	double *sample_best;
//	double **sample_time_series;
//	double **sample_dist;
	int sample_idx;
	double *best_v;
	int *best_m_l;
	int * sample_series_index;
	best_v = (double*) malloc(1*sizeof(double));
	best_m_l = (int*) malloc(2*sizeof(int));
	
	B = get_block_size(time_begin, time_end);
	L = (time_end-time_begin)/B + 1;
//    cout<<"B, L"<<B<<','<<L<<endl;
	sample_best=(double*) malloc(R*sizeof(double));

	sample_series_index = (int *) malloc((time_end - time_begin) * sizeof(int));
    srand(123);
	for(r=0;r<R;r++)
	{	
//		cout<<"r: "<<r<<endl;
		for(l=0;l<L;l++)
		{
			sample_idx =  time_begin + rand()%(time_end-time_begin - B + 1); 
			for(b=0;b<B;b++)
			{	
				if(B*l+b<time_end-time_begin)
				{
				    sample_series_index[B*l+b] = sample_idx+b;
				}
				else
				{
					break;
				}		
			}
			if(B*l+b>=time_end-time_begin)
			{
				break; 
			}
		}
		powell(best_m_l, best_v, 0, time_end-time_begin, sample_series_index);
		sample_best[r]=best_v[0];		
	}

	int sort_num= int(R * pt);
	double *max_sort;
	double threshold;
	max_sort=(double*) malloc(sort_num*sizeof(double));
	for(i=0;i<sort_num;i++)
	{
	    max_sort[i] = 0;
	}


	for(r=0;r<R;r++)
	{
//	    cout<<"sample_best: "<<sample_best[r]<<endl;
		for(i=0;i<sort_num;i++)
		{
			if(sample_best[r] > max_sort[i])
			{
			    if(sort_num > 1)
			    {
			        for(j=sort_num-2;j>=i;j--)
				    {
					    max_sort[j + 1] = max_sort[j];
				    }
				    max_sort[i] = sample_best[r];
			    }
			    else
			    {
			        max_sort[0] = sample_best[r];
			    }
				break;
			}

		}
//		cout<<"max_sort: "<<max_sort[sort_num - 1]<<endl;
	}


	threshold=max_sort[sort_num - 1];
//    cout<<"threshold end 0"<<endl;
    free(best_v);
	free(best_m_l);
    free(max_sort);
    free(sample_best);
    free(sample_series_index);
//    cout<<"threshold end 1"<<endl;
//	free_matrix(sample_time_series, time_end-time_begin, time_series_d);
//    cout<<"threshold end 2"<<endl;
//	free_matrix(sample_dist, time_end - time_begin, time_end - time_begin);
//    cout<<"threshold end 3"<<endl;
	return threshold;
}

void cp::cp_origin()
{
	change_point_set[0] = 0;
	change_point_set[1] = time_series_len;
	change_point_num[0] = 0;
	
	double max_V_each_segment[time_series_len];
	int max_m_l_each_segment[time_series_len];	
	int new_segment_index[2];	
	new_segment_index[0] = 0;
	
	int new_segment_num=1;
	
	int s, i, m, l;
	int s_num, time_begin, time_end;
	double v;

	int* series_index;
	series_index = (int *) malloc(time_series_len * sizeof(int));
	for(i=0;i<time_series_len;i++)
	{
	    series_index[i] = i;
	}
	while(1)
	{ 
		for(s_num=0;s_num<new_segment_num;s_num++)
		{
			s = new_segment_index[s_num];		
			time_begin=change_point_set[s];
			time_end=change_point_set[s+1];
//			cout<<"time_begin: "<<time_begin<<endl;
//			cout<<"time_end: "<<time_end<<endl;

			if(time_end - time_begin <= 2 * min_segment_size)
			{
				for(i=time_series_len;i>s;i--)
				{
					max_V_each_segment[i] = max_V_each_segment[i - 1];
					max_m_l_each_segment[i] = max_m_l_each_segment[i - 1];		
				}
				max_V_each_segment[i] = 0;
				max_m_l_each_segment[i] = -1;
			}
			else
			{
				double best=0;
				int best_m=0;
				for(m = time_begin + min_segment_size; m<time_end - min_segment_size; m++)
				{
					for(l = m + 1; l<time_end; l++)
					{
						v = calculate_V(series_index, time_begin, m, l);
						if(v > best)
						{
							best = v;
							best_m = m;
						}
					}
				}
				for(i=time_series_len;i>s;i--)
				{
					max_V_each_segment[i] = max_V_each_segment[i - 1];
					max_m_l_each_segment[i] = max_m_l_each_segment[i - 1];		
				}
				max_V_each_segment[i] = best;
				max_m_l_each_segment[i] = best_m;
			}
			
		}
		
		double max_max_V=-1;
		int max_max_V_index=0;
		int max_time_begin=0;
		int max_time_end=0;
		int max_m;
		
		for(i=0;i<change_point_num[0]+1;i++)
		{
			if(max_V_each_segment[i]>max_max_V)
			{
				max_time_begin=change_point_set[i];
				max_time_end=change_point_set[i+1];
				max_max_V = max_V_each_segment[i];
				max_max_V_index=i;
			}		
		}
//		cout<<"max_max_V: " <<max_max_V<<endl;
		double threshold = cal_threshold_cp(max_time_begin, max_time_end);
//		cout<<"threshold: " <<threshold<<endl;

		if(max_max_V > threshold)
		{
			max_m = max_m_l_each_segment[max_max_V_index];
			for(i=change_point_num[0]+2;i>max_max_V_index+1;i--)
			{
				change_point_set[i]=change_point_set[i-1];
				
			}
			change_point_set[i]=max_m;
			change_point_num[0] += 1;
			
			for(i=max_max_V_index;i<change_point_num[0] -1;i++)
			{
				max_V_each_segment[i]=max_V_each_segment[i+1];
				max_m_l_each_segment[i]=max_m_l_each_segment[i+1];			
			}		
			new_segment_index[0] = max_max_V_index;
			new_segment_index[1] = max_max_V_index+1;
			new_segment_num = 2;		
		}
		else
		{
            free(series_index);
            free_matrix(time_series,time_series_len,time_series_d);
            free_matrix(dist,time_series_len,time_series_len);
			return;
		}		
	}	

}

void cp::cp_gss()
{
	int s_num, time_begin, time_end, s, i, m;
	change_point_set[0] = 0;
	change_point_set[1] = time_series_len;
	change_point_num[0] = 0;
	
	double max_V_each_segment[time_series_len];
	int max_m_l_each_segment[time_series_len];
	int new_segment_index[2];	
	new_segment_index[0] = 0;
	int new_segment_num=1;
	double best;

	int *best_m_l;
	double *this_m_best_v;
	int *this_m_best_m_l;
	best_m_l =(int*) malloc(2*sizeof(int));
	this_m_best_v =(double*) malloc(1*sizeof(double));
	this_m_best_m_l =(int*) malloc(2*sizeof(int));
	
	int *p, *u;
	p =(int*) malloc(2*sizeof(int));
	u =(int*) malloc(2*sizeof(int));
	double threshold;

	int* series_index;
	series_index = (int *) malloc(time_series_len * sizeof(int));
	for(i=0;i<time_series_len;i++)
	{
	    series_index[i] = i;
	}

    double max_max_V = -1;
    int max_max_V_index;
    int max_time_begin;
    int max_time_end;
    int max_m;
	
	while(1)
	{ 
		for(s_num=0;s_num<new_segment_num;s_num++)
		{
			s = new_segment_index[s_num];		
			time_begin=change_point_set[s];
			time_end=change_point_set[s+1];
//			cout<<"time_begin: "<<time_begin<<endl;
//			cout<<"time_end: "<<time_end<<endl;
//
			if(time_end - time_begin <= 2 * min_segment_size)
			{
				for(i=time_series_len;i>s;i--)
				{
					max_V_each_segment[i] = max_V_each_segment[i - 1];
					max_m_l_each_segment[i] = max_m_l_each_segment[i - 1];		
				}
				max_V_each_segment[i] = 0;
				max_m_l_each_segment[i] = -1;
			}

			else
			{	
				best = 0.;
				best_m_l[0]=0;
				best_m_l[1]=0;
				for(m = time_begin + min_segment_size; m<time_end - min_segment_size; m++)
				{
                    p[0] = m;
					p[1] = m+1;
					u[0] = 0;
					u[1] = 1;

					golden_section_search(this_m_best_m_l, this_m_best_v, p, u, time_begin, time_end, series_index, 2);

					if(this_m_best_v[0] > best)
					{
						best = this_m_best_v[0];
						best_m_l[0] = this_m_best_m_l[0];
						best_m_l[1] = this_m_best_m_l[1];						
					}
//					if(m == time_begin + min_segment_size)
//                    {
//                        cout<<7;
//                    }
				}

//				cout<<2;

				for(i=time_series_len;i>s;i--)
				{
					max_V_each_segment[i] = max_V_each_segment[i - 1];
					max_m_l_each_segment[i] = max_m_l_each_segment[i - 1];		
				}
				max_V_each_segment[i] = best;
				max_m_l_each_segment[i] = best_m_l[0];
			}
		}

		max_max_V=-1;
        max_max_V_index=0;
        max_time_begin=0;
        max_time_end=0;
        max_m=0;

		for(i=0;i<change_point_num[0]+1;i++)
		{
			if(max_V_each_segment[i]>max_max_V)
			{
				max_time_begin=change_point_set[i];
				max_time_end=change_point_set[i+1];
				max_max_V = max_V_each_segment[i];
				max_max_V_index=i;
			}		
		}

		threshold = cal_threshold_cp_gss(max_time_begin, max_time_end);

		if(max_max_V > threshold)
		{
//            cout<<2;
			max_m = max_m_l_each_segment[max_max_V_index];
			for(i=change_point_num[0]+2;i>max_max_V_index+1;i--)
			{
				change_point_set[i]=change_point_set[i-1];
				
			}
			change_point_set[i]=max_m;
			change_point_num[0] += 1;
			
			for(i=max_max_V_index;i<change_point_num[0] -1;i++)
			{
				max_V_each_segment[i]=max_V_each_segment[i+1];
				max_m_l_each_segment[i]=max_m_l_each_segment[i+1];			
			}		
			new_segment_index[0] = max_max_V_index;
			new_segment_index[1] = max_max_V_index+1;
			new_segment_num = 2;
//			cout<<3;
		}
		else
		{
            free(this_m_best_m_l);
            free(best_m_l);
            free(this_m_best_v);
            free(p);
            free(u);
            free(series_index);
            free_matrix(time_series,time_series_len,time_series_d);
            free_matrix(dist,time_series_len,time_series_len);
			return;
		}		
	}	
}

void cp::cp_powell()
{	
	change_point_set[0] = 0;
	change_point_set[1] = time_series_len;
	change_point_num[0]=0;
	
	double max_V_each_segment[time_series_len];
	int max_m_l_each_segment[time_series_len];	
	int new_segment_index[2];	
	new_segment_index[0] = 0;
	
	int new_segment_num=1;
	int s, i;

	int* series_index;
	series_index = (int *) malloc(time_series_len * sizeof(int));
	for(i=0;i<time_series_len;i++)
	{
	    series_index[i] = i;
	}
	

	int s_num, time_begin, time_end;
	
	double *best_v;
	int *best_m_l;
	
	best_v = (double*) malloc(1*sizeof(double));
	best_m_l = (int*) malloc(2*sizeof(int));
	
	while(1)
	{ 
		for(s_num=0;s_num<new_segment_num;s_num++)
		{
			s = new_segment_index[s_num];		
			time_begin=change_point_set[s];
			time_end=change_point_set[s+1];
//			cout<<"time_begin: "<<time_begin<<endl;
//			cout<<"time_end: "<<time_end<<endl;
			
			if(time_end - time_begin <= 2 * min_segment_size)
			{
				for(i=time_series_len;i>s;i--)
				{
					max_V_each_segment[i] = max_V_each_segment[i - 1];
					max_m_l_each_segment[i] = max_m_l_each_segment[i - 1];		
				}
				max_V_each_segment[i] = 0;
				max_m_l_each_segment[i] = -1;
			}
			else
			{
//				cout<<"powell"<<endl;
				powell(best_m_l, best_v, time_begin, time_end, series_index);
//				cout<<"1"<<endl;
				for(i=time_series_len;i>s;i--)
				{
					max_V_each_segment[i] = max_V_each_segment[i - 1];
					max_m_l_each_segment[i] = max_m_l_each_segment[i - 1];		
				}
				max_V_each_segment[i] = best_v[0];
				max_m_l_each_segment[i] = best_m_l[0];
			}
			
		}
		
		double max_max_V=-1;
		int max_max_V_index=0;
		int max_time_begin=0;
		int max_time_end=0;
		int max_m;
		for(i=0;i<change_point_num[0]+1;i++)
		{
			if(max_V_each_segment[i]>max_max_V)
			{
				max_time_begin=change_point_set[i];
				max_time_end=change_point_set[i+1];
				max_max_V = max_V_each_segment[i];
				max_max_V_index=i;
			}		
		}
//		cout<<"max_max_V: "<<max_max_V<<endl;
//		cout<<"max_time_begin: "<<max_time_begin<<endl;
//		cout<<"max_time_end: "<<max_time_end<<endl;
		double threshold = cal_threshold_cp_powell(max_time_begin, max_time_end);
//		cout<<"threshold: "<<threshold<<endl;

		if(max_max_V > threshold)
		{

			max_m = max_m_l_each_segment[max_max_V_index];
			for(i=change_point_num[0]+2;i>max_max_V_index+1;i--)
			{
				change_point_set[i]=change_point_set[i-1];
				
			}
			change_point_set[i]=max_m;
			change_point_num[0] += 1;
			
			for(i=max_max_V_index;i<change_point_num[0];i++)
			{
				max_V_each_segment[i]=max_V_each_segment[i+1];
				max_m_l_each_segment[i]=max_m_l_each_segment[i+1];			
			}		
			new_segment_index[0] = max_max_V_index;
			new_segment_index[1] = max_max_V_index+1;
			new_segment_num = 2;		
		}
		else
		{
            free(best_v);
            free(best_m_l);
            free(series_index);
            free_matrix(time_series,time_series_len,time_series_d);
            free_matrix(dist,time_series_len,time_series_len);
			return;
		}		
	}	
}

void run(double *my_time_series, int my_time_series_len, int my_time_series_d, int my_min_segment_size, double my_pt, bool my_is_Euclidean_distance, double *my_dist, int my_block_size, int * my_change_point_set, int * my_change_point_num, char method)
{
    cp cp_1(my_time_series, my_time_series_len, my_time_series_d, my_min_segment_size, my_pt, my_is_Euclidean_distance, my_dist, my_block_size, my_change_point_set, my_change_point_num);
    if(method == 'o')
    {
        cp_1.cp_origin();
    }
    else if(method == 'g')
    {
        cp_1.cp_gss();
    }
    else
    {
        cp_1.cp_powell();
    }
}
 

//int main()
//{
//	int i, j;
//	double *my_time_series;
//	int ts_len=80;
//	int ts_d=3;
//
//	my_time_series=(double*) malloc((ts_len*ts_d)*sizeof(double));
//
//	for(i=0;i<40;i++)
//	{
//		for(j=0;j<ts_d;j++)
//		{
//			my_time_series[i*ts_d + j]=1;
//		}
//	}
//
//	my_time_series[20] = 2.0;
//
//	for(i=40;i<80;i++)
//	{
//		for(j=0;j<ts_d;j++)
//		{
//			my_time_series[i*ts_d+j]=5;
//		}
//	}
//	cp cp_1(my_time_series, ts_len, ts_d, 3, 0.05, '1');
//
//	cout<<"init finished"<<endl;
//
//	int *change_point_set, *change_point_num;
//	change_point_set = (int*) malloc(ts_len*sizeof(int));
//	change_point_num = (int*) malloc(1*sizeof(int));
////	cp_1.cp_origin(change_point_set, change_point_num);
////	cp_1.cp_gss(change_point_set, change_point_num);
//	cp_1.cp_powell(change_point_set, change_point_num);
//}
