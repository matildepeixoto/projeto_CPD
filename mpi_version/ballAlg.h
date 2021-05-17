#ifndef FUNC_H
#define FUNC_H

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p)) 
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1)

struct node {
  long id;
  double *center;
  long leaf;
  double radius;
  int *procsList;
  int sizeList;
  struct node *nextL; 
  struct node *nextR;
};

struct node *newNode;

//#define PRINT_TIME2 1


double **create_array_pts(int n_dims, long np);
double **get_points(int argc, char *argv[], int *n_dims, long *np);
double get_distance(int n_dims, double *x1, double *x2);
void get_points_ab_mpi(double **pts, long *set, int n_dims, long n_points, long *a, long *b, int id, int p, MPI_Comm comm);
void get_points_ab(double **pts, long *set, int n_dims, long n_points, long *a, long *b);
void orthogonal_projection_mpi(double **pts, long *set, double **po, int n_dims, long n_points, long a, long b, int id, int p, MPI_Comm comm);
void orthogonal_projection(double **pts, long *set, double **po, int n_dims, long n_points, long a, long b);
void calc_median(double **pts, long *set, int n_dims, long i, long a, long b, double *median);
static int comp(const void *p1, const void *p2);
void find_median(double **pts, long *set, double **po, int n_dims, long n_points, long a, long b, double *median);
double get_radius(double **pts, long *set, long n_points, int n_dims, double *median);
void create_sets_LR(long *set, double **po, int n_dims, long n_points, double *median, long *l, long *r);
struct node* build_tree_mpi(double **pts, long *set, double **po, int n_dims, long n_points, int gid, int id, int gp, int p, MPI_Comm comm, int *procsList);
struct node* build_tree(double **pts, long *set, double **po, int n_dims, long n_points, int *procsList);
void freepointers(long n, double** pointer);

#endif