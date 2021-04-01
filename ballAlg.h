#ifndef FUNC_H
#define FUNC_H
/*
struct node {
  double *center;
  double radius;
  struct node *next;
};

struct node *newNode;

struct Tree {
};
*/
double **create_array_pts(int n_dims, long np);
double **get_points(int argc, char *argv[], int *n_dims, long *np);
double get_distance(int n_dims, double *x1, double *x2);
void get_points_ab(double **pts, int n_dims, long n_points, int *a, int *b);
double **orthogonal_projection(double **pts, int n_dims, long n_points, int a, int b);
static int comp(const void *p1, const void *p2);
double *find_median(double **po, int n_dims, long n_points);
double get_radius(double **pts, int n_dims, int a, int b, double *median);
void create_sets_LR(double **pts, double **set_L, double **set_R, double **po, int n_dims, int n_points, double *median);

#endif