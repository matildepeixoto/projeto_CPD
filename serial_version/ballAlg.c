#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ballAlg.h"

#define RANGE 10

long node_id = 0;

struct node* createNode(int n_dims) {

    struct node* newNode = malloc(sizeof(struct node));
    newNode->center = (double*) malloc(n_dims * sizeof(double));
    newNode->nextL = NULL;
    newNode->nextR = NULL;

    return newNode;
}

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
    
    for(int i = 0; i < np; i++){
        p_arr[i] = &_p_arr[i * n_dims];
    }

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
        dist = dist + (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }

    dist = sqrt(dist);

    return dist;
}

void get_points_ab(double **pts, long *set, int n_dims, long n_points, long *a, long *b)
{
    long i, aux = 0;
    double dist = 0, dist_aux;

    aux = set[0];

    for(i = 1; i < n_points; i++)
    {
        dist_aux = get_distance(n_dims, pts[aux], pts[set[i]]);
        if(dist_aux > dist)
        {
            dist = dist_aux;
            *a = i;
        }
    }
    
    dist = 0;
    aux = set[*a];

    for(i = 0; i < n_points; i++)
    {
        dist_aux = get_distance(n_dims, pts[aux], pts[set[i]]);
        if(dist_aux > dist)
        {
            dist = dist_aux;
            *b = i;
        }
    }
}

void orthogonal_projection(double **pts, long *set, double **po, int n_dims, long n_points, long a, long b)
{
    double den = 0, num = 0, inn_prod = 0, aux = 0;
    double *aux2 = (double *) malloc(n_dims * sizeof(double));
    long index_a = set[a], index_b = set[b], index_i = 0;

    for(int j = 0; j < n_dims; j++)
    {
        aux = pts[index_b][j] - pts[index_a][j];
        den += aux * aux;
        aux2[j] = aux;
    }
    
    for(long i = 0; i < n_points; i++)
    {
        if(i != a && i != b)
        {   
            index_i = set[i];
            for(int j = 0; j < n_dims; j++)
            {
                num += (pts[index_i][j] - pts[index_a][j]) * aux2[j];
            }
            inn_prod = (num/den);
            po[i][0] = (inn_prod * (pts[index_b][0] - pts[index_a][0])) + pts[index_a][0];
            num = 0;
        }
        else if(i == a)
            po[i][0] = pts[index_a][0];
        else
            po[i][0] = pts[index_b][0];
        po[i][1] = set[i];
    }

    free(aux2);
}

void calc_median(double **pts, long *set, int n_dims, long i, long a, long b, double *median)
{
    
    double den = 0, num = 0, inn_prod = 0, aux = 0;
    double *aux2 = (double *) malloc(n_dims * sizeof(double));
    long index_a = set[a], index_b = set[b];

    for(int j = 0; j < n_dims; j++)
    {
        aux = pts[index_b][j] - pts[index_a][j];
        den += aux * aux;
        aux2[j] = aux; 
    }

    if(i != index_a && i != index_b)
    {
        for(int j = 0; j < n_dims; j++)
        {
            num += (pts[i][j] - pts[index_a][j]) * aux2[j];
        }

        inn_prod = (num/den);
        
        for(int j = 0; j < n_dims; j++)
        {
            median[j] = (inn_prod * aux2[j]) + pts[index_a][j];
        }   
        num = 0;
    }
    else if(i == index_a)
    {
        for(int j = 0; j < n_dims; j++)
        {
            median[j] = pts[index_a][j];
        }
    }
    else
    {
        for(int j = 0; j < n_dims; j++)
        {
            median[j] = pts[index_b][j];
        }
    }

    free(aux2);
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

void find_median(double **pts, long *set, double **po, int n_dims, long n_points, long a, long b, double *median)
{
    int index = 0, idx1 = 0, idx2 = 0;
    
    qsort(po, n_points, sizeof(po[0]), comp);

    if(n_points%2 == 1)
    {        
        index = n_points/2;
        calc_median(pts, set, n_dims, po[index][1], a, b, median);
    }
    else
    {
        double *median_aux = (double*)malloc(n_dims * sizeof(double));
        idx1 = n_points/2;
        idx2 = idx1 - 1;
        calc_median(pts, set, n_dims, po[idx1][1], a, b, median);
        calc_median(pts, set, n_dims, po[idx2][1], a, b, median_aux);
        for (int i = 0; i < n_dims; i++)
            median[i] = (median[i] + median_aux[i])/2;
        free(median_aux);
    }
}

double get_radius(double **pts, long *set, long n_points, int n_dims, double *median)
{
    double dist_aux = 0, radius = 0;
    
    for(long i = 0; i < n_points; i++)
    {
        dist_aux = get_distance(n_dims, pts[set[i]], median);
        if(dist_aux > radius)
            radius = dist_aux;
    }

    return radius;
}

void create_sets_LR(long *set, double **po, int n_dims, long n_points, double *median, long *l, long *r)
{    
    long l_aux = 0, r_aux = n_points/2;

    for(long i = 0; i < n_points; i++)
    {   
        if(po[i][0] < median[0])
        {
            set[l_aux] = po[i][1];
            l_aux++;
        }
        else
        {
            set[r_aux] = po[i][1];
            r_aux++;
        }
    }   
    *l = l_aux;
    *r = r_aux - n_points/2;
}

struct node* build_tree(double **pts, long *set, double **po, int n_dims, long n_points)
{   
    double radius = 0;
    struct node* root;

    root = createNode(n_dims);

    root->id = node_id;
    node_id++; 
    
    if(n_points > 1)
    {
        long a, b, l = 0, r = 0;
        double dist;

        get_points_ab(pts, set, n_dims, n_points, &a, &b);
        
        orthogonal_projection(pts, set, po, n_dims, n_points, a, b);

        find_median(pts, set, po, n_dims, n_points, a, b, root->center);

        root->radius = get_radius(pts, set, n_points, n_dims, root->center);

        create_sets_LR(set, po, n_dims, n_points, root->center, &l, &r);      
        
        root->nextL = build_tree(pts, set, po, n_dims, l);
        root->nextR = build_tree(pts, &set[l], &po[l], n_dims, r);       
    }
    else
    {
        root->radius = 0;
        for(int i = 0; i < n_dims; i++)
            root->center[i] = pts[set[0]][i];
    }
    
    return root;    
}

void print_tree(struct node* Tree, int n_dims)
{
    struct node* tempL = Tree;
    struct node* tempR = Tree;
    
    printf("%ld", tempL->id);
    if(tempL->nextL == NULL)
        printf(" -1");
    else
        printf("  %ld", (tempL->nextL)->id);
    if(tempL->nextR == NULL)
        printf(" -1");
    else
        printf("  %ld", (tempL->nextR)->id);
    printf("  %.6lf", tempL->radius);
    for(int i = 0; i < n_dims; i++)
        printf("  %.6lf", tempL->center[i]);
    printf("\n");
    tempL = tempL->nextL;
    if (tempL != NULL)
        print_tree(tempL, n_dims);
    else
    {
        free(Tree->center);
        free(Tree);
        return;
    }
    tempR = tempR->nextR;
    if (tempR != NULL)
        print_tree(tempR, n_dims);
    else
    {
        free(Tree->center);
        free(Tree);
        return;
    }
    free(Tree->center);
    free(Tree);
}

void freepointers(long n, double** pointer){
    for (long i = 0; i < n; i++)
    {
        free(pointer[i]);
    }
    free(pointer);
}

int main(int argc, char *argv[])
{
    double exec_time;
    double **pts;
    int n_dims;
    long n_points;
    struct node* root;

    exec_time = -omp_get_wtime();
    pts = get_points(argc, argv, &n_dims, &n_points);
    long* set = (long*)malloc(n_points * sizeof(long)); 
    double** po = (double**)malloc((n_points) * sizeof(double*));
    for (long i = 0; i < n_points; i++)
    {
        set[i] = i;
        po[i] = (double*)malloc(2 * sizeof(double));
    }
    root = build_tree(pts, set, po, n_dims, n_points);
    exec_time += omp_get_wtime();
    free(pts[0]);
    free(pts);
    freepointers(n_points, po);
    fprintf(stderr, "%.10lf\n", exec_time);
    printf("%d %ld\n", n_dims, node_id);
    //print_tree(root, n_dims);
}