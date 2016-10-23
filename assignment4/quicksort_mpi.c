#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

typedef struct node_ {
    int key;
    int value;
} node_t;
typedef struct heap_ {
    node_t *arr;
    int size;
    int heap_size;
} heap_t;
#define PARENT(i) i/2
#define LEFT(i) 2*i
#define RIGHT(i) 2*i+1
#define MAX_VALUE INT_MAX
#define TRUE 1
#define FALSE 0
void swap(heap_t heap, int i, int j)
{
    node_t temp;
    temp = heap.arr[i];
    heap.arr[i] = heap.arr[j];
    heap.arr[j] = temp;
    return;
}
void min_heapify(heap_t heap, int i)
{
    int l = LEFT(i);
    int r = RIGHT(i);
    int smallest = i;
    if (l <= heap.heap_size && heap.arr[l].key < heap.arr[smallest].key) {
	smallest = l;
    }
    if (r <= heap.heap_size && heap.arr[r].key < heap.arr[smallest].key) {
	smallest = r;
    }
    if (smallest != i) {
	swap (heap, i, smallest);
	min_heapify(heap, smallest);
    }
}
void build_min_heap(heap_t heap)
{
    int i = 0;
    for (i = heap.heap_size/2; i>=1; i--) {
	min_heapify(heap, i);
    }
}
int heap_min(heap_t heap)
{
    return heap.arr[1].key;
}
node_t extract_min_heap(heap_t *heap)
{
    node_t temp;
    if (heap->heap_size == 0) {
	temp.key = MAX_VALUE;
	temp.value = MAX_VALUE;
	return temp;
    }
    temp = heap->arr[1];
    heap->arr[1] = heap->arr[heap->heap_size];
    heap->heap_size -= 1;
    min_heapify(*heap, 1);
    return temp;
}
void heap_decrease_key(heap_t heap, int i, int key)
{
    if (heap.arr[i].key < key) {
	printf("Value is less than existing value");
	return;
    }
    heap.arr[i].key = key;
    while (i >= 1 && PARENT(i) >= 1
	   && heap.arr[i].key < heap.arr[PARENT(i)].key) {
	swap(heap, i, PARENT(i));
	i = PARENT(i);
    }
    return;
}
void heap_insert(heap_t *heap, node_t key)
{
    heap->heap_size += 1;
    if (heap->heap_size >= heap->size) {
	printf("Over flow\n");
	return;
    }
    heap->arr[heap->heap_size].key = MAX_VALUE;
    heap->arr[heap->heap_size].value = key.value;
    heap_decrease_key(*heap, heap->heap_size, key.key);
}
int is_underflow(node_t node)
{
    if (node.key == MAX_VALUE && node.value == MAX_VALUE) {
	return TRUE;
    }
    return FALSE;
}
int merge_sorted_lists(int *result, heap_t *sorted, int k, int *n)
{
    int tot = 0;
    heap_t min_heap;
    node_t elem, dummy;
    int i = 0;
    min_heap.arr = malloc((k+1)*sizeof(node_t));
    min_heap.size = k+1;
    min_heap.heap_size = k;
    for (i= 1; i<=k ; i++) {
	tot += sorted[i-1].heap_size;
	min_heap.arr[i] = extract_min_heap(&sorted[i-1]);
    }
    
    build_min_heap(min_heap);
    for (i = 0; i<tot; i++) {
	elem = extract_min_heap(&min_heap);
	result[i] = elem.key;
	dummy = extract_min_heap(&sorted[elem.value]);
	if (is_underflow(dummy) == TRUE) {
	    continue;
	}
	heap_insert(&min_heap, dummy);
    }
    *n = tot;
    return 0;
}
void create_min_heap(heap_t *sorted, int *arr, int n, int k)
{
    int i = 0;
    sorted[k].arr = malloc((n+1)*sizeof(node_t));
    sorted[k].size = n+1;
    sorted[k].heap_size = n;
    sorted[k].arr[0].key = MAX_VALUE;
    sorted[k].arr[0].value = k;
    for (i = 1; i<=n; i++) {
	sorted[k].arr[i].key = arr[i-1];
	sorted[k].arr[i].value = k;
    }
}
void alloc_array(int **arr, int N);
void fill_array(int **arr, int N, int num_procs);
int cmpfunc(const void *a, const void *b);
int main (int argc, char *argv[])
{
    const int N = 64;
    int rank, num_procs;
    int *arr = NULL;
    int *s_disp = NULL;
    int *s_count = NULL;
    int *pivot_buf = NULL;
    int *s_res_start = NULL;
    int *s_res_len = NULL;
    int *r_buf = NULL;
    int *r_len = NULL;
    int *r_start = NULL;
    int *temp_buf = NULL;
    
    int root = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    alloc_array(&arr, N);
    alloc_array(&s_disp, num_procs);
    alloc_array(&s_count, num_procs);
    alloc_array(&pivot_buf, num_procs*num_procs);
    alloc_array(&temp_buf, num_procs*num_procs);

    for (int i = 0; i<num_procs; i++) {
	s_count[i] = N/num_procs;
	s_disp[i] = i*s_count[i];
    }
    if (rank == root) {
	fill_array(&arr, N, num_procs);
    }
    for (int i = 1; (rank == root) && i<N ;i++) {
	/* if (arr[i] < arr[i-1]) { */
	/*     printf("Not Sorted\n"); */
	/*     break; */
	/* } */
	printf("%d\n", arr[i]);
    }
    if (rank == root) {
	MPI_Scatterv(arr, s_count, s_disp, MPI_INT, MPI_IN_PLACE,
		     s_count[rank], MPI_INT, root, MPI_COMM_WORLD);
    } else  {
	MPI_Scatterv(arr, s_count, s_disp, MPI_INT, arr,
		     s_count[rank], MPI_INT, root, MPI_COMM_WORLD);
    }
    qsort(arr, s_count[rank], sizeof(int), cmpfunc);
    /* if (rank == 1) { */
    /* 	for (int i = 0; i<s_count[rank]; i++) { */
    /* 	    printf("Rank: %d Data: %d\n", rank, arr[i]); */
    /* 	} */
    /* } */
    for(int i=0;i<num_procs;i++) {
	pivot_buf[i]= arr[i*s_count[rank]/num_procs];
    }
    /* if (rank == 1) { */
    /* 	for (int i = 0; i<num_procs; i++) { */
    /* 	    printf("Rank: %d Data: %d\n", rank, pivot_buf[i]); */
    /* 	} */
    /* } */
    if(rank == root) {
	MPI_Gather(MPI_IN_PLACE, num_procs, MPI_INT,
		   pivot_buf, num_procs, MPI_INT, root, MPI_COMM_WORLD);
    } else {
	MPI_Gather(pivot_buf, num_procs, MPI_INT,
		   pivot_buf, num_procs, MPI_INT, root, MPI_COMM_WORLD);
    }

    if (rank == root) {
	heap_t *hp = malloc(num_procs*sizeof(heap_t));
	int l2;
	for (int i = 0; i < num_procs; i++) {
	    create_min_heap(hp, &pivot_buf[i*num_procs], num_procs, i);
	}
	merge_sorted_lists(temp_buf, hp, num_procs, &l2);
	for(int i=0; i<num_procs-1; i++) {
	    pivot_buf[i] = temp_buf[(i+1)*num_procs];
	}
    }
    if (rank == root) {
    	for (int i = 0; i<num_procs*num_procs; i++) {
    	    printf("Rank: %d Data: %d\n", rank, pivot_buf[i]);
    	}
    }
    MPI_Bcast(pivot_buf, num_procs-1, MPI_INT, root, MPI_COMM_WORLD);
    alloc_array(&s_res_len, num_procs);
    alloc_array(&s_res_start, num_procs);
    int k = 0;
    for (int i = 0; i<num_procs-1; i++) {
	s_res_len[i] = k;
	s_res_start[i] = 0;

	while ((k < s_count[rank]) && (arr[k] <= pivot_buf[i])) {
	    s_res_len[i]++;
	    k++;
	}
    }
    s_res_start[num_procs - 1] = k;
    s_res_len[num_procs - 1] = s_count[rank] - k;

    alloc_array(&r_buf, N);
    alloc_array(&r_len, num_procs);
    alloc_array(&r_start, num_procs);

    for (int i = 0; i < num_procs; i++) {
	MPI_Gather(&s_res_len[i], 1, MPI_INT, 
		   r_len, 1, MPI_INT, i, MPI_COMM_WORLD);
	if (rank == i) {
	    r_start[0] = 0;
	    for (int j = 0; j < num_procs; j++) {
		r_start[j] = r_start[j] + r_len[j-1];
	    }
	}
	MPI_Gatherv(&arr[s_res_start[i]],
		    s_res_len[i], MPI_INT,
		    r_buf, r_len, r_start, MPI_INT, i, MPI_COMM_WORLD);
    }
    heap_t *sorted = malloc(num_procs*sizeof(heap_t));
    int l;
    for (int i = 0; i < num_procs; i++) {
	create_min_heap(sorted, r_buf+r_start[i], r_len[i], i);
    }
    merge_sorted_lists(arr, sorted, num_procs, &l);

    int send_len = r_start[num_procs-1] + r_len[num_procs-1];
    int send_lens[num_procs];
    int send_starts[num_procs];

    MPI_Gather(&send_len, 1, MPI_INT,
	       send_lens, 1 ,MPI_INT, root, MPI_COMM_WORLD);

    if (rank == root) {
	send_starts[0] = 0;
	for(int i = 1; i<num_procs; i++) {
	    send_starts[i] = send_starts[i-1] + send_lens[i-1];
	}	
    }

    int *res = NULL;
    alloc_array(&res, N);
    MPI_Gatherv(arr, send_len, MPI_INT,
		res, send_lens, send_starts, MPI_INT, root, MPI_COMM_WORLD);

    MPI_Finalize();
    if (rank == root) {
	printf("\n\n");
    }
    for (int i = 0; (rank == root) && i<N ;i++) {
	/* if (res[i] < res[i-1]) { */
	/*     printf("Not Sorted\n"); */
	/*     return 0; */
	/* } */
	printf("%d\n", res[i]);
    }
    if (rank == root) {
	printf("Sorted\n");
    }
    return 0;
}

void alloc_array(int **arr, int N)
{
    *arr = malloc(N*sizeof(int));
}
void fill_array(int **arr, int N, int num_procs)
{
    int *arr_l = *arr;
    int i = 0;

    for (i = 0; i<N; i++)
	arr_l[i] = rand()%1000;
}
int cmpfunc(const void *a, const void *b)
{
    /*Copied from 
      https://www.tutorialspoint.com/c_standard_library/c_function_qsort.htm*/
    return (*(int*)a - *(int*)b);
}
