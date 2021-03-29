#ifndef FUNC_H
#define FUNC_H

double **create_array_pts(int n_dims, long np);
double **get_points(int argc, char *argv[], int *n_dims, long *np);
double get_distance(int n_dims, double *x1, double *x2);
void get_points_ab(double **pts, int n_dims, long n_points, int *a, int *b);
double **orthogonal_projection(double **pts, int n_dims, long n_points, int a, int b);
#endif