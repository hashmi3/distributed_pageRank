#ifndef UTILIT_H
#define UTILIT_H
//#define __STDC_FORMAT_MACROS 1
#include <inttypes.h>
//#include <cmath>



uintV countTriangles(uintV *array1, uintE len1, uintV *array2, uintE len2,
                     uintV u, uintV v) {

  uintE i = 0, j = 0; // indexes for array1 and array2
  uintV count = 0;

  if (u == v)
    return count;

  while ((i < len1) && (j < len2)) {
    if (array1[i] == array2[j]) {
      if ((array1[i] != u) && (array1[i] != v)) {
        count++;
      }
      i++;
      j++;
    } else if (array1[i] < array2[j]) {
      i++;
    } else {
      j++;
    }
  }
  return count;
}

void get_thread_load(int *thrd_load, size_t workers, size_t n){

    int rem = n%workers;
    
    for(uint i=0; i< workers; i++){
        thrd_load[i] = n/workers;
        if(rem > 0){
               thrd_load[i] += 1;
               rem--;
        }
    }        
}



void get_load_4_equal_edge(Graph &g, int* thrd_load, int* thrd_edge, int n_workers, uintV n ){
    uintE tot_edg = 0;
    uintE temp = 0;
  
    uint * vert_edg_cnt = (uint *) calloc(n,sizeof(uint));

    for(uint i=0; i< n; i++){
        temp = g.vertices_[i].getOutDegree();
        vert_edg_cnt[i] += temp;
        tot_edg += temp;
    } 
    //printf("Total edges: %" PRId32 "\n edge Div: ", tot_edg );
    get_thread_load(thrd_edge, n_workers, tot_edg);   //get equal edge division for threads
  
    int j =0;
    int ptr= 0;
    temp =0;
    for(uint i=0; i<n; i++){
        temp += vert_edg_cnt[i];
        if (temp >= thrd_edge[j]){
            //printf("temp: %lu diff: %lu, ", temp, temp-thrd_edge[j]);
            thrd_load[j] = i - ptr;
            ptr = i;
            temp = 0;
            j++;
        }
    }
    thrd_load[j] = n-ptr;
    
    free(vert_edg_cnt);
}


void get_load_4_equal_edge_pr(Graph &g, uintV* thrd_load, uintV* thrd_edge, int n_workers, uintV n /*, double time_prt*/ ){
    //timer t_prt;
    uintE tot_edg = 0;
    uintE temp = 0;
  
    uint * vert_edg_cnt = (uint *) calloc(n,sizeof(uint));

    //t_prt.start();
    for(uint i=0; i< n; i++){
        temp = g.vertices_[i].getOutDegree();
        vert_edg_cnt[i] += temp;
        tot_edg += temp;
    } 
    //printf("Total edges: %" PRId32 "\n edge Div: ", tot_edg );
    get_thread_load(thrd_edge, n_workers, tot_edg);   //get equal edge division for threads
  
    int j =0;
    int ptr= 0;
    temp =0;
    for(uint i=0; i<n; i++){
        temp += vert_edg_cnt[i];
        if (temp >= thrd_edge[j]){
            //printf("temp: %lu diff: %lu, ", temp, temp-thrd_edge[j]);
            thrd_load[j] = i - ptr;
            ptr = i;
            temp = 0;
            j++;
        }
    }
    thrd_load[j] = n-ptr;
    //time_prt = t_prt.stop();

    free(vert_edg_cnt);
}


void calcu_1(Graph &g, int id, uintV v_start, uintV v_end, PageRankType *pr_curr, PageRankType *pr_holdr, PageRankType *pr_next, int max_iters /*, int * pr_visit*/, PageRankType * thrd_sum, uintV n, uint workers, uintV *split_ptr_begin, uintV *thrd_load, PageRankType &tot_sum ){

    timer t1, t_b1, t_b2;
    double time_taken = 0.0;
    double time_b1 = 0.0;
    double time_b2 = 0.0;
    uintE out_degree;
    uintV v;
    uintE edge_count =0;
    uint cnt_vert = 0;
    int world_rank;
    PageRankType local_sum = 0;

    
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size( MPI_COMM_WORLD, &world_size);
    
    int sendcounts[ world_size ];
    int displs[ world_size ];
    
    t1.start();
    
    for (int iter = 0; iter < max_iters; iter++) {
    // for each vertex 'u', process all its outNeighbors 'v'
    for (uintV u = v_start; u < v_end; u++) {
      out_degree = g.vertices_[u].getOutDegree();
      edge_count += out_degree;
      for (uintE i = 0; i < out_degree; i++) {
        v = g.vertices_[u].getOutNeighbor(i);
        pr_next[v] += (pr_curr[u] / out_degree);               
      }
    }
    t_b1.start();
    
    MPI_Reduce( pr_next, pr_holdr, n, PAGERANK_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD );
    if( id == 0){
      for (uint i=0; i< n; i++){
        pr_next[i]= pr_holdr[i];
      }
    }
    MPI_Scatterv( pr_next, thrd_load, split_ptr_begin, PAGERANK_MPI_TYPE, pr_holdr , thrd_load[id], PAGERANK_MPI_TYPE, 0, MPI_COMM_WORLD);

    time_b1 += t_b1.stop();
    
    //printf("Pass Barr ID: %d \n",id);
    for (uintV v = v_start; v < v_end; v++) {
      pr_next[v] = PAGE_RANK(pr_holdr[v - split_ptr_begin[id] ]);           //this eliminates above for loop

      // reset pr_curr for the next iteration
      pr_curr[v] = pr_next[v];
      cnt_vert++;
    }
    for (uintV i=0; i< n; i++){
      pr_next[i] = 0.0;
    }
    
    }
    time_taken = t1.stop();
    printf("%d, %d, %lf\n",id, edge_count, time_b1);

    // Each process calculates partial sum for its vertx range 
    for(uintV i=v_start; i<v_end; i++ ){
      local_sum += pr_curr[i];
    }
    
    //Non root processes send local_sum to root
    if (id != 0)
      MPI_Send( &local_sum , 1, PAGERANK_MPI_TYPE, 0, 0, MPI_COMM_WORLD );
    //Root accumulates local_sum and add its own
    PageRankType tmp_count = 0;
    if(id == 0){
      for( uint i=1; i< workers; i++ ){
        MPI_Recv( &tmp_count, 1, PAGERANK_MPI_TYPE , i, 0, MPI_COMM_WORLD , MPI_STATUS_IGNORE);
        tot_sum += tmp_count;
        
      }
      tot_sum +=local_sum;    
    }
    
}


void calcu_2(Graph &g, int id, uintV v_start, uintV v_end, PageRankType *pr_curr, PageRankType *pr_holdr, PageRankType *pr_next, int max_iters /*, int * pr_visit*/, PageRankType * thrd_sum, uintV n, uint workers, uintV *split_ptr_begin, uintV *thrd_load, PageRankType &tot_sum ){

    timer t1, t_b1, t_b2;
    double time_taken = 0.0;
    double time_b1 = 0.0;
    double time_b2 = 0.0;
    uintE out_degree;
    uintV v;
    uintE edge_count =0;
    uint cnt_vert = 0;
    int world_rank;
    PageRankType local_sum = 0;

    
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size( MPI_COMM_WORLD, &world_size);
    
    int sendcounts[ world_size ];
    int displs[ world_size ];
    
    t1.start();
    for (int iter = 0; iter < max_iters; iter++) {
    // for each vertex 'u', process all its outNeighbors 'v'
    for (uintV u = v_start; u < v_end; u++) {
      out_degree = g.vertices_[u].getOutDegree();
      edge_count += out_degree;
      for (uintE i = 0; i < out_degree; i++) {
        v = g.vertices_[u].getOutNeighbor(i);
        pr_next[v] += (pr_curr[u] / out_degree);               
      }
    }
    t_b1.start();

    for(int i=0; i< world_size; i++){
      MPI_Reduce( pr_next + split_ptr_begin[i]  , pr_holdr , thrd_load[i] , PAGERANK_MPI_TYPE, MPI_SUM, i, MPI_COMM_WORLD );
    }

    time_b1 += t_b1.stop();
    
    for (uintV v = v_start; v < v_end; v++) {
      pr_next[v] = PAGE_RANK(pr_holdr[v - split_ptr_begin[id] ]);           

      // reset pr_curr for the next iteration
      pr_curr[v] = pr_next[v];
      //pr_next[v] = 0.0;
      cnt_vert++;
    }
    for (uintV i=0; i< n; i++){
      pr_next[i] = 0.0;
    }
    }
    time_taken = t1.stop();
    printf("%d, %d, %lf\n",id, edge_count, time_b1);


    // Each process calculates partial sum for its vertx range 
    for(uintV i=v_start; i<v_end; i++ ){
      local_sum += pr_curr[i];
    }
    //Non root processes send local_sum to root
    if (id != 0)
      MPI_Send( &local_sum , 1, PAGERANK_MPI_TYPE, 0, 0, MPI_COMM_WORLD );
    //Root accumulates local_sum and add its own
    PageRankType tmp_count = 0;
    if(id == 0){
      for( uint i=1; i< workers; i++ ){
        MPI_Recv( &tmp_count, 1, PAGERANK_MPI_TYPE , i, 0, MPI_COMM_WORLD , MPI_STATUS_IGNORE);
        tot_sum += tmp_count;
      }
      tot_sum +=local_sum;    
    }
}




#endif

