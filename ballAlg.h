#ifndef FUNC_H
#define FUNC_H

struct node {
  long id;
  double *center;
  double radius;
  struct node *nextL; 
  struct node *nextR;
};

struct node *newNode;

//#define PRINT_TIME2 1


double **create_array_pts(int n_dims, long np);
double **get_points(int argc, char *argv[], int *n_dims, long *np);
double get_distance(int n_dims, double *x1, double *x2);
void get_points_ab(double **pts, long *set, int n_dims, long n_points, long *a, long *b);
void orthogonal_projection(double **pts, long *set, double **po, int n_dims, long n_points, long a, long b);
double *calc_median(double **pts, long *set, int n_dims, long i, long a, long b);
static int comp(const void *p1, const void *p2);
double *find_median(double **pts, long *set, double **po, int n_dims, long n_points, long a, long b);
double get_radius(double **pts, long *set, int n_points, int n_dims, double *median);
void create_sets_LR(long *set_L, long *set_R, double **po, int n_dims, long n_points, double *median, long *l, long *r);
void freepointers(long n, double** pointer);

#endif