#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ballAlg.h"

#define RANGE 10
/*
struct node* createNode(int n_dims, double *median, double radius) {

    struct node* newNode = malloc(sizeof(struct node));
    newNode->center = (double*) malloc(n_dims * sizeof(double));

    for(int i = 0; i < n_dims; i++)
        newNode->center[i] = median[i];
    newNode->radius = radius;
    newNode->next = NULL;

    return newNode;
}

// Create a graph of V vertices
struct Tree* create_tree() {

  struct GrapTreeh* tree = malloc(sizeof(struct Tree));

  tree->num_V = 0;

  // Create vertical array of nodes (size MAX_NODES)
  graph->a_list = (struct node **) malloc(MAX_NODES * sizeof(struct node *));
  graph->visited =(int*) malloc(MAX_NODES * sizeof(int*));
  graph->tier1 = (int*) malloc(MAX_NODES * sizeof(int*));

  //initializes the vectors
  for(int i = 0; i < MAX_NODES ; i++ )
  {
    graph->visited[i] = 0;
    graph->a_list[i] = NULL;
    graph->tier1[i] = 0;
  }

  return graph;
}
*/
double **create_array_pts(int n_dims, long np)
{
    double *_p_arr;
    double **p_arr;
    
    _p_arr = (double *) malloc(n_dims * np * sizeof(double));
    p_arr = (double **) malloc(np * sizeof(double *));
    
    if((_p_arr == NULL) || (p_arr == NULL)){
        printf("Error allocating array of points, exiting.\n");
        exit(4);
    }
    
    for(int i = 0; i < np; i++)
        p_arr[i] = &_p_arr[i * n_dims];
    
    return p_arr;
}

double **get_points(int argc, char *argv[], int *n_dims, long *np)
{
    double **pt_arr;
    unsigned seed;
    long i;
    int j;

    if(argc != 4){
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }

    *n_dims = atoi(argv[1]);

    if(*n_dims < 2){
        printf("Illegal number of dimensions (%d), must be above 1.\n", *n_dims);
        exit(2);
    }

    *np = atol(argv[2]);

    if(*np < 1){
        printf("Illegal number of points (%ld), must be above 0.\n", *np);
        exit(3);
    }
    
    seed = atoi(argv[3]);
    srandom(seed);
    pt_arr = create_array_pts(*n_dims, *np);

    for(i = 0; i < *np; i++){
        for(j = 0; j < *n_dims; j++)
            pt_arr[i][j] = RANGE * ((double) random()) / RAND_MAX;
    }
    
    return pt_arr;
}

double get_distance(int n_dims, double *x1, double *x2)
{
    int i = 0;
    double dist = 0;

    for(i = 0; i < n_dims; i++)
    {
        dist = dist + pow((x1[i] - x2[i]), 2);
    }

    dist = sqrt(dist);

    return dist;
}

void get_points_ab(double **pts, int n_dims, long n_points, int *a, int *b)
{
    long i;
    double dist = 0, dist_aux;
    
    for(i = 1; i < n_points; i++)
    {
        dist_aux = get_distance(n_dims, pts[0], pts[i]);
        if(dist_aux > dist)
        {
            dist = dist_aux;
            *a = i;
        }
    }
    
    dist = 0;

    for(i = 0; i < n_points; i++)
    {
        dist_aux = get_distance(n_dims, pts[*a], pts[i]);
        if(dist_aux > dist)
        {
            dist = dist_aux;
            *b = i;
        }
    }
}

double **orthogonal_projection(double **pts, int n_dims, long n_points, int a, int b)
{
    double den, num, inn_prod;
    double** po = (double**)malloc((n_points) * sizeof(double*));

    for (int i = 0; i < n_points; i++)
    {
        po[i] = (double*)malloc((n_dims+1) * sizeof(double));
    }
    
    for(int i = 0; i < n_points; i++)
    {
        if(i != a && i != b)
        {
            for(int j = 0; j < n_dims; j++)
            {
                num += (pts[i][j] - pts[a][j]) * (pts[b][j] - pts[a][j]);
                den += pow((pts[b][j] - pts[a][j]), 2);
            }
            inn_prod = (num/den);
            for(int j = 0; j < n_dims; j++)
            {
                po[i][j] = (inn_prod * (pts[b][j] - pts[a][j])) + pts[a][j];
            }
            num = 0;
            den = 0;
        }
        else if(i == a)
        {
            for(int j = 0; j < n_dims; j++)
            {
                po[a][j] = pts[a][j];
            }
        }
        else
        {
            for(int j = 0; j < n_dims; j++)
            {
                po[b][j] = pts[b][j];
            }
        }
        po[i][n_dims] = i;
    }

    return po;
}

static int comp(const void *p1, const void *p2) {
    const double (*x1) = *(const double**)p1;
    const double (*x2) = *(const double**)p2;

    double diff = x1[0] - x2[0];
    if (diff > 0)
        return 1;
    else if (diff < 0)
        return -1;
    else
        return 0;
}

double *find_median(double **po, int n_dims, long n_points)
{
    int index = 0, idx1 = 0, idx2 = 0;
    double *median = (double*)malloc(n_dims * sizeof(double));

    qsort(po, n_points, sizeof(po[0]), comp);

    if(n_points%2 == 1)
    {
        index = div(n_points, 2).quot;
        printf("index: %d\n", index);
        for(int i = 0; i < n_dims; i++)
            median[i] = po[index][i];
    }
    else
    {
        idx1 = div(n_points, 2).quot;
        idx2 = idx1 + 1;
        for(int i = 0; i < n_dims; i++)
        {
            median[i] = (po[idx1][i] + po[idx2][i])/2;
        }
    }

    return median;
}

double get_radius(double **pts, int n_dims, int a, int b, double *median)
{
    double dist_a = 0, dist_b = 0, radius = 0;

    dist_a = get_distance(n_dims, pts[a], median);
    dist_b = get_distance(n_dims, pts[b], median);
    if(dist_a > dist_b)
        radius = dist_a;
    else
        radius = dist_b;

    return radius;
}

void create_sets_LR(double **pts, double **set_L, double **set_R, double **po, int n_dims, int n_points, double *median)
{
    int l = 0, r = 0;
    
    for(int i = 0; i < n_points; i++)
    {
        if(po[i][0] < median[0])
        {
            for(int j = 0; j < n_dims; j++)
                set_L[l][j] = pts[((int)po[i][n_dims])][j];
            l++;
        }
        else
        {
            for(int j = 0; j < n_dims; j++)
                set_R[r][j] = pts[((int)po[i][n_dims])][j];
            r++;
        }
    }   
}

int main(int argc, char *argv[])
{
    double exec_time, radius;
    double *median;
    double **pts, **po;
    int n_dims;
    long n_points;
    int a, b;
    double dist;

    exec_time = -omp_get_wtime();
    pts = get_points(argc, argv, &n_dims, &n_points);
    int points_set = n_points/2 + 1;
    double** set_L = (double**)malloc(points_set * sizeof(double*));
    double** set_R = (double**)malloc(points_set * sizeof(double*));
    
    for (int i = 0; i < n_points; i++)
    {
        set_L[i] = (double*)malloc((n_dims) * sizeof(double));
        set_R[i] = (double*)malloc((n_dims) * sizeof(double));
    }
    get_points_ab(pts, n_dims, n_points, &a, &b);
    printf("Pontos:\n");
    for(int i = 0; i < n_points; i++)
        printf("%.2lf %.2lf\n", pts[i][0], pts[i][1]);
    printf("a: %.10lf, %.10lf;\nb: %.10lf, %.10lf\n", pts[a][0], pts[a][1], pts[b][0], pts[b][1]);
    po = orthogonal_projection(pts, n_dims, n_points, a, b);
    printf("Projeção ortogonal:\n");
    for(int i = 0; i < n_points; i++)
        printf("%.2lf %.2lf\n", po[i][0], po[i][1]);
    median = find_median(po, n_dims, n_points);
    printf("Projeção ortogonal ordenada:\n");
    for(int i = 0; i < n_points; i++)
        printf("%.2lf %.2lf\n", po[i][0], po[i][1]);
    printf("Median:\n");
    for(int i = 0; i < n_dims; i++)
        printf("%.2lf\n", median[i]);
    radius = get_radius(pts, n_dims, a, b, median);
    printf("Radius: %.2lf\n", radius);
    create_sets_LR(pts, set_L, set_R, po, n_dims, n_points, median);
    printf("chegueei\n");
    printf("Set L:\n");
    for(int i = 0; i < 2; i++)
        printf("%.2lf %.2lf\n", set_L[i][0], set_L[i][1]); 
    printf("Set R:\n");
    for(int i = 0; i < 3; i++)
        printf("%.2lf %.2lf\n", set_R[i][0], set_R[i][1]); 
    //root = build_tree();
    exec_time += omp_get_wtime();
    fprintf(stderr, "%.10lf\n", exec_time);
    //dump_tree(root); // to the stdout!
}