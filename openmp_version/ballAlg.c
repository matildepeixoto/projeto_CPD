#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ballAlg.h"

#define RANGE 10

//OMP_NUM_THREADS=8;

long node_id = -1;

struct node* createNode(int n_dims, double *median, double radius, long id) {

    struct node* newNode = malloc(sizeof(struct node));
    newNode->center = (double*) malloc(n_dims * sizeof(double));

    newNode->id = id;
    for(int i = 0; i < n_dims; i++)
        newNode->center[i] = median[i];
    newNode->radius = radius;
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

    #pragma omp parallel for if(n_points > 10000) private(dist_aux) 
    for(i = 1; i < n_points; i++)
    {
        dist_aux = get_distance(n_dims, pts[aux], pts[set[i]]);
        #pragma omp critical (find_a)
        if(dist_aux > dist)
        {
            dist = dist_aux;
            *a = i;
        }
    }
    
    dist = 0;
    aux = set[*a];

    #pragma omp parallel for if(n_points > 10000) private(dist_aux) 
    for(i = 0; i < n_points; i++)
    {
        dist_aux = get_distance(n_dims, pts[aux], pts[set[i]]);
        #pragma omp critical (find_b)
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
    
    #pragma omp parallel for if(n_points > 10000) private(index_i, inn_prod) firstprivate(num)
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
            po[a][0] = pts[index_a][0];
        else
            po[b][0] = pts[index_b][0];
        po[i][1] = set[i];
    }
    free(aux2);
}

double *calc_median(double **pts, long *set, int n_dims, long i, long a, long b)
{
    double den = 0, num = 0, inn_prod = 0, aux = 0;
    double *median = (double*)malloc(n_dims * sizeof(double));
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

    return median;
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

double *find_median(double **pts, long *set, double **po, int n_dims, long n_points, long a, long b)
{
    int index = 0, idx1 = 0, idx2 = 0;
    double *median, *median_aux;

    qsort(po, n_points, sizeof(po[0]), comp);

    if(n_points%2 == 1)
    {        
        index = n_points/2;
        median = calc_median(pts, set, n_dims, po[index][1], a, b);
    }
    else
    {
        idx1 = n_points/2;
        idx2 = idx1 - 1;
        median = calc_median(pts, set, n_dims, po[idx1][1], a, b);
        median_aux = calc_median(pts, set, n_dims, po[idx2][1], a, b);
        for (int i = 0; i < n_dims; i++)
            median[i] = (median[i] + median_aux[i])/2;
    }

    return median;
}

double get_radius(double **pts, long *set, int n_points, int n_dims, double *median)
{
    double dist_aux = 0, radius = 0;
    
    #pragma omp parallel for if(n_points > 10000) private(dist_aux) //schedule(dynamic)
    for(long i = 0; i < n_points; i++)
    {
        dist_aux = get_distance(n_dims, pts[set[i]], median);
        #pragma omp critical (find_radius)
        if(dist_aux > radius)
            radius = dist_aux;
    }

    return radius;
}

void create_sets_LR(long *set_L, long *set_R, double **po, int n_dims, long n_points, double *median, long *l, long *r)
{    
    long l_aux = 0, r_aux = 0, aux = 0;

    for(long i = 0; i < n_points; i++)
    {   
        if(po[i][0] < median[0])
        {
            set_L[l_aux] = po[i][1];
            l_aux++;
        }
        else
        {
            set_R[r_aux] = po[i][1];
            r_aux++;            
        }
    }   
    *l = l_aux;
    *r = r_aux;
}

struct node* build_tree(double **pts, long *set, int n_dims, long n_points)
{   
    long id;
    double *median = pts[set[0]];
    double radius = 0;
    struct node* root;

    #pragma omp critical
    {
        node_id++; 
        id = node_id;
    }
    
    if(n_points > 1)
    {
        long a, b, l = 0, r = 0;
        double dist;
        

        get_points_ab(pts, set, n_dims, n_points, &a, &b);
        
        orthogonal_projection(pts, set, po, n_dims, n_points, a, b);
        
        median = find_median(pts, set, po, n_dims, n_points, a, b);
        
        #pragma omp parallel sections if(n_points > 1000)
        {
            #pragma omp section
            radius = get_radius(pts, set, n_points, n_dims, median);
            #pragma omp section
            create_sets_LR(set_L, set_R, po, n_dims, n_points, median, &l, &r);
        }
        
        root = createNode(n_dims, median, radius, id);
        
        #pragma omp parallel 
        {
            #pragma omp sections nowait
            {
                #pragma omp section
                root->nextL = build_tree(pts, set_L, n_dims, l);
                #pragma omp section
                root->nextR = build_tree(pts, set_R, n_dims, r);
            }
        }
        free(set_L);
        free(set_R);
        freepointers(n_points, po);
        free(median);
    }
    else
    {
        root = createNode(n_dims, median, radius, id); 
        
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
    /*long* set_L = (long*)malloc(n_points * sizeof(long));
    long* set_R = (long*)malloc(n_points * sizeof(long));*/
    double** po = (double**)malloc((n_points) * sizeof(double*));   
    for (long i = 0; i < n_points; i++)
    {
        po[i] = (double*)malloc(2 * sizeof(double));
        set[i] = i;
    }
    root = build_tree(pts, set, n_dims, n_points);
    exec_time += omp_get_wtime();
    free(pts[0]);
    free(pts);
    fprintf(stderr, "%.10lf\n", exec_time);
    printf("%d %ld\n", n_dims, node_id);
    print_tree(root, n_dims);
}