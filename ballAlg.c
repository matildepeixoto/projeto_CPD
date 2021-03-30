#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ballAlg.h"

#define RANGE 10

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
    double** po = (double**)malloc(n_points * sizeof(double*));

    for (int i = 0; i < n_points; i++)
    {
        po[i] = (double*)malloc(n_dims * sizeof(double));
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

int main(int argc, char *argv[])
{
    double exec_time;
    double *median;
    double **pts, **po;
    int n_dims;
    long n_points;
    int a, b;
    double dist;

    exec_time = -omp_get_wtime();
    pts = get_points(argc, argv, &n_dims, &n_points);
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
    //root = build_tree();
    exec_time += omp_get_wtime();
    fprintf(stderr, "%.10lf\n", exec_time);
    //dump_tree(root); // to the stdout!
}