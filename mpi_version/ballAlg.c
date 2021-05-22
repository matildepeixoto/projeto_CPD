#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include <unistd.h>
#include <string.h>
#include "ballAlg.h"

#define RANGE 10

long node_id = 0;

struct node* createNode(int n_dims, int p) {

    struct node* newNode = malloc(sizeof(struct node));
    newNode->center = (double*) malloc(n_dims * sizeof(double));
    newNode->procsList = (int*) malloc(p * sizeof(int));
    newNode->leaf = -1;
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
        sleep(2);
        exit(1);
    }

    *n_dims = atoi(argv[1]);

    if(*n_dims < 2){
        printf("Illegal number of dimensions (%d), must be above 1.\n", *n_dims);
        sleep(2);
        exit(2);
    }

    *np = atol(argv[2]);

    if(*np < 1){
        printf("Illegal number of points (%ld), must be above 0.\n", *np);
        sleep(2);
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

    return dist;
}

void get_points_ab_mpi(double **pts, long *set, int n_dims, long n_points, long *a, long *b, int id, int p, MPI_Comm comm)
{
    long i = 0, aux = 0, gi = 0;
    double dist = 0, dist_aux, gdist = 0;
    MPI_Status status;

    aux = set[0];

    *a = 0;
    *b = 1;

    long size = BLOCK_SIZE(id, p, n_points);
    for(i = 0; i < size; i++)
    {
        gi = i + BLOCK_LOW(id, p, n_points);
        dist_aux = get_distance(n_dims, pts[aux], pts[set[gi]]);
        if(dist_aux > dist)
        {
            dist = dist_aux;
            *a = gi;
        }
    }

    MPI_Allreduce(&dist, &gdist, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(gdist == dist)
    {
        for(i = 0; i < p; i++)
        {
            if(i != id)
                MPI_Send(a, 1, MPI_LONG, i, 1, comm);
        }
        
    }
    else
    {
        MPI_Recv(a, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, comm, &status);
    }
    
    dist = 0;
    aux = set[*a];

    for(i = 0; i < size; i++)
    {
        gi = i + BLOCK_LOW(id, p, n_points);
        dist_aux = get_distance(n_dims, pts[aux], pts[set[gi]]);
        if(dist_aux > dist)
        {
            dist = dist_aux;
            *b = gi;
        }
    }

    MPI_Allreduce(&dist, &gdist, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(gdist == dist)
    {
        for(i = 0; i < p; i++)
        {
            if(i != id)
                MPI_Send(b, 1, MPI_LONG, i, 1, comm);
        }       
    }
    else
    {
        MPI_Recv(b, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, comm, &status);
    }
}

void get_points_ab(double **pts, long *set, int n_dims, long n_points, long *a, long *b)
{
    long i, aux = 0, index;
    double dist = 0, dist_aux;

    aux = set[0];

    *a = 0;
    *b = 1;

    for(i = 1; i < n_points; i++)
    {
        index = set[i];
        dist_aux = get_distance(n_dims, pts[aux], pts[index]);
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
        index = set[i];
        dist_aux = get_distance(n_dims, pts[aux], pts[index]);
        if(dist_aux > dist)
        {
            dist = dist_aux;
            *b = i;
        }
    }
}
void orthogonal_projection_mpi(double **pts, long *set, double **po, double **aux_array, int n_dims, long n_points, long a, long b, int id, int p,  MPI_Comm comm)
{
    double den = 0, num = 0, inn_prod = 0, aux = 0;
    double *aux2 = (double *) malloc(n_dims * sizeof(double));
    long index_a = set[a], index_b = set[b], index_i = 0, gi = 0;
    long size = BLOCK_SIZE(id, p, n_points), pos_init = BLOCK_LOW(id, p, n_points);

    int * recvcounts = (int *) malloc (p * sizeof(int));
    int * displs = (int *) malloc (p * sizeof(int));
    for (int i = 0; i < p; i++)
    {
        recvcounts[i] = 2 * BLOCK_SIZE(i, p, n_points);
        displs[i] = 2 * BLOCK_LOW(i, p, n_points);
    }


    for(int j = 0; j < n_dims; j++)
    {
        aux = pts[index_b][j] - pts[index_a][j];
        den += aux * aux;
        aux2[j] = aux;
    }

    for(long i = 0; i < size; i++)
    {
        gi = i + pos_init;

        if(gi != a && gi != b)
        {   
            index_i = set[gi];
            for(int j = 0; j < n_dims; j++)
            {
                num += (pts[index_i][j] - pts[index_a][j]) * aux2[j];
            }
            inn_prod = (num/den);
            aux_array[i][0] = (inn_prod * (pts[index_b][0] - pts[index_a][0])) + pts[index_a][0];
            num = 0;
        }
        else if(gi == a)
            aux_array[i][0] = pts[index_a][0];
        else
            aux_array[i][0] = pts[index_b][0];
        aux_array[i][1] = set[gi];
    }
        

    MPI_Allgatherv(&(aux_array[0][0]), 2 * size, MPI_DOUBLE, &(po[0][0]), recvcounts, displs, MPI_DOUBLE, comm);

    free(aux2);
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
            po[a][0] = pts[index_a][0];
        else
            po[b][0] = pts[index_b][0];
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

/*procedure quicksort(array, left, right)
{

    if(right > left)
    {

    }
        select a pivot index (e.g. pivotIndex := left)
        pivotNewIndex := partition(array, left, right, pivotIndex)
        quicksort(array, left, pivotNewIndex - 1)
        quicksort(array, pivotNewIndex + 1, right)
}

function partition(array, left, right, pivotIndex)
{
    pivotValue := array[pivotIndex]
    swap array[pivotIndex] and array[right] // Move pivot to end
    storeIndex := left
    for i from left to right - 1
        if array[i] <= pivotValue
            swap array[i] and array[storeIndex]
            storeIndex := storeIndex + 1
    // Move pivot to its final place
    swap array[storeIndex] and array[right]

    return storeIndex;
}*/

void swap(double* a, double* b)
{
    double t[2];

    t[0] = a[0];
    t[1] = a[1];
    
    a[0] = b[0];
    a[1] = b[1];

    b[0] = t[0];
    b[1] = t[1];
}

long partition (double **arr, long low, long high)
{
    double pivot = arr[high][0]; // pivot
    long i = (low - 1); // Index of smaller element and indicates the right position of pivot found so far
 
    for (long j = low; j <= high - 1; j++)
    {
        // If current element is smaller than the pivot
        if (arr[j][0] < pivot)
        {
            i++; // increment index of smaller element
            swap(&arr[i][0], &arr[j][0]);
        }
    }
    swap(&arr[i + 1][0], &arr[high][0]);
    return (i + 1);
}

void quickSort(double **arr, long low, long high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
        at right place */
        long pi = partition(arr, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

void find_median(double **pts, long *set, double **po, int n_dims, long n_points, long a, long b, double *median)
{
    int index = 0, idx1 = 0, idx2 = 0;
    
    //qsort(po, n_points, sizeof(po[0]), comp);
    quickSort(po, 0, n_points-1);

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
    radius = sqrt(radius);
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

struct node* build_tree_mpi(double **pts, long *set, double **po, double **aux_array, int n_dims, long n_points, int gid, int id, int gp, int p, MPI_Comm comm, int *procsList)
{   
    int new_id, new_p, id_bcast;
    double radius = 0;
    struct node* root;
    
    
    //MPI_Status status;

    root = createNode(n_dims, p);
    
    for(int i = 0; i < p; i++)
        root->procsList[i] = procsList[i];
    root->sizeList = p;
    //root->id = node_id;
    if(gid == root->procsList[0])
    {
        node_id++;   
        /*for(int i = 0; i < gp; i++)
        {
            if(i != gid)
                MPI_Send(&node_id, 1, MPI_LONG, i, 4, MPI_COMM_WORLD);
        } */    
    }
    /*else
    {
        MPI_Recv(&node_id, 1, MPI_LONG, MPI_ANY_SOURCE, 4, MPI_COMM_WORLD, &status);
    }*/
    
    if(n_points > 1)
    {
        long a, b, l = 0, r = 0;
        double dist;
        
        get_points_ab_mpi(pts, set, n_dims, n_points, &a, &b, id, p, comm);

        

        //orthogonal_projection_mpi(pts, set, po, aux_array, n_dims, n_points, a, b, id, p, comm);


        /*if ( id == 0 )
        {*/

            orthogonal_projection(pts, set, po, n_dims, n_points, a, b);


            find_median(pts, set, po, n_dims, n_points, a, b, root->center);

/* 
            for (int i = 0; i < n_points; i++)
            {
                printf("%d &(po[%d][0]): %p &(po[%d][1]): %p \n", gid, i, &(po[i][0]), i, &(po[i][1]));
            }
             for (int i = 0; i < n_points; i++)
            {
                printf("%d &(test[%d][0]): %p &(test[%d][1]): %p \n", gid, i, &(test[i][0]), i, &(test[i][1]));
            } */

            root->radius = get_radius(pts, set, n_points, n_dims, root->center);

        
            create_sets_LR(set, po, n_dims, n_points, root->center, &l, &r); 
            
            
        //}
        if(p > 1)
        {
            MPI_Comm new_comm; 
            int color = id/(p/2);
            MPI_Comm_split(comm, color, id, &new_comm);
            MPI_Comm_rank(new_comm, &new_id);
            MPI_Comm_size(new_comm, &new_p);
            if(id < p/2){
                //printf("left id: %d %d \n", id, new_id);
                root->nextL = build_tree_mpi(pts, set, po, aux_array, n_dims, l, gid, new_id, gp, new_p, new_comm, procsList);
            }
            else{
                //printf("right id: %d %d \n", id, new_id);
                root->nextR = build_tree_mpi(pts, &set[l], &po[l], &aux_array[l], n_dims, r, gid, new_id, gp, new_p, new_comm, &procsList[new_p]);  
            }
        }
        else
        {
            #pragma omp task shared(root) //final(2^level > num_t)
            root->nextL = build_tree(pts, set, po, n_dims, l, procsList);
            #pragma omp task shared(root) //final(2^level > num_t)
            root->nextR = build_tree(pts, &set[l], &po[l], n_dims, r, procsList);
        }
         
    }
    else
    {
        root->radius = 0;
        root->leaf = set[0];
    }
    
    return root;    
}

struct node* build_tree(double **pts, long *set, double **po, int n_dims, long n_points, int *procsList)
{   
    struct node* root;
    
    root = createNode(n_dims, 1);

    #pragma omp critical
    {
        root->procsList[0] = procsList[0];
        root->sizeList = 1;
        node_id++; 
        //MPI_Bcast(&node_id, 1, MPI_LONG, procsList[0], MPI_COMM_WORLD);
        root->id = node_id;
    }
    
    if(n_points > 1)
    {
        long a, b, l = 0, r = 0;
        double dist;
              
        get_points_ab(pts, set, n_dims, n_points, &a, &b);
        
        orthogonal_projection(pts, set, po, n_dims, n_points, a, b);
        
        find_median(pts, set, po, n_dims, n_points, a, b, root->center);

        root->radius = get_radius(pts, set, n_points, n_dims, root->center);
        
        create_sets_LR(set, po, n_dims, n_points, root->center, &l, &r);

        #pragma omp task shared(root) //final(2^level > num_t)
        root->nextL = build_tree(pts, set, po, n_dims, l, procsList);
        #pragma omp task shared(root) //final(2^level > num_t)
        root->nextR = build_tree(pts, &set[l], &po[l], n_dims, r, procsList); 
    }
    else
    {
        root->radius = 0;
        root->leaf = set[0];
    }  

    return root;    
}

void print_tree(struct node* Tree, int n_dims, double** pts, int id)
{
    int check = 0;
    long id_parent = 0;
    char aux[1000] = "", str[1000] = "";
    struct node* tempL = Tree;
    struct node* tempR = Tree;

    id_parent = tempL->id;
    if(id == tempL->procsList[0])
    {
        if(tempL->nextL == NULL){
            sprintf(str, "%ld", tempL->id);
            strcpy(aux, " -1 -1 ");
            strcat(str, aux);
            sprintf(aux, "%.6lf", tempL->radius);
            strcat(str, aux);
            //printf("%ld -1 -1 %.6lf", tempL->id, tempL->radius);
            for(int i = 0; i < n_dims; i++){
                sprintf(aux, "  %.6lf", pts[tempL->leaf][i]);
                strcat(str, aux);
               // printf("  %.6lf", pts[tempL->leaf][i]);
            }
        }           
            
        else{
            sprintf(aux, "%ld", tempL->id);
            strcpy(str, aux);
            sprintf(aux, " %ld", 2*id_parent + 1);
            strcat(str, aux);
            sprintf(aux, " %ld", 2*id_parent + 2);
            strcat(str, aux);
            sprintf(aux, " %.6lf", tempL->radius);
            strcat(str, aux);
            /* printf("%ld %ld %ld %.6lf", tempL->id, 2*id_parent + 1, 2*id_parent + 2, tempL->radius); */
            for(int i = 0; i < n_dims; i++){
                sprintf(aux, "  %.6lf", tempL->center[i]);
                strcat(str, aux);
                /* printf("  %.6lf", tempL->center[i]); */
            }
        }
        strcat(str, "\n");
        printf("%s", str);
        //printf("\n");
    }
    tempL = tempL->nextL;
    if(tempL != NULL)
    {
        tempL->id = 2*id_parent + 1;
        print_tree(tempL, n_dims, pts, id);
    }
    tempR = tempR->nextR; 
    if(tempR != NULL)
    {
        tempR->id = 2*id_parent + 2;
        print_tree(tempR, n_dims, pts, id);
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
    int n_dims, me, nprocs;
    long n_points, total_nodes;
    struct node* root;

    MPI_Init(&argc, &argv);
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    exec_time = - MPI_Wtime();

    int *procsList = (int*)malloc(nprocs * sizeof(int)); ;
    for(int i = 0; i < nprocs; i++)
        procsList[i] = i;

    pts = get_points(argc, argv, &n_dims, &n_points);
    long* set = (long*)malloc(n_points * sizeof(long)); 

    double** po = (double **)malloc(n_points * sizeof(double*));
    //create_array_pts(2, n_points);
    double *data = (double *)malloc(2*n_points*sizeof(double));
    double** aux_array = create_array_pts(2, n_points);

    //double **po = (double **)malloc(n_points * sizeof(double *));
    for (long i = 0; i < n_points; i++)
    {
        set[i] = i;
        po[i] = &data[i * 2];
        
    }
    root = build_tree_mpi(pts, set, po, aux_array, n_dims, n_points, me, me, nprocs, nprocs, MPI_COMM_WORLD, procsList);

    exec_time += MPI_Wtime();
    
    MPI_Reduce(&node_id, &total_nodes, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(me == 0)
    {
        fprintf(stderr, "%.1lf\n", exec_time);
        printf("%d %ld\n", n_dims, total_nodes);
  
    } 
    MPI_Barrier(MPI_COMM_WORLD);
    root->id = 0;
    print_tree(root, n_dims, pts, me); 
    
    //nao funciona com free po[0]??????????
    free(data);
    free(po);
    free(pts[0]);
    free(pts);

    MPI_Finalize();
}