#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <thread>
#include <mutex>
#include <stdexcept>
#define __STDC_FORMAT_MACROS 1
#include <inttypes.h>
#include "utilit.h"
#include <mpi.h>
#include <assert.h>


 #ifdef USE_INT
 #define INIT_PAGE_RANK 100000
 #define EPSILON 1000
 #define PAGE_RANK(x) (15000 + (5 * x) / 6)
 #define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
 #define PAGERANK_MPI_TYPE MPI_LONG
 #define PR_FMT "%ld"
 typedef int64_t PageRankType;
 #else
 #define INIT_PAGE_RANK 1.0
 #define EPSILON 0.01
 #define DAMPING 0.85
 #define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
 #define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
 #define PAGERANK_MPI_TYPE MPI_FLOAT
 #define PR_FMT "%f"
 typedef float PageRankType;
 #endif





int pageRankParallel_2(Graph &g, int workers, int world_rank, int max_iters, uintV * thrd_load, uintV *thrd_edge ) {

  uintV n = g.n_;

  PageRankType *pr_curr = new PageRankType[n];
  PageRankType *pr_next = new PageRankType[n];
  PageRankType *pr_holdr = new PageRankType[n];
  
  timer t1;
  double time_taken = 0.0;
  double time_prt = 0.0;
  PageRankType tot_sum = 0;
  PageRankType local_sum = 0;

    timer  t_b1;
    double time_b1 = 0.0;
    uintE out_degree;
    uintV v;
    uintE edge_count =0;
    uintV v_start, v_end;

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
    pr_holdr[i] = 0.0;
  }

  //int world_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size = workers;
  int id = world_rank;
  //MPI_Comm_size(MPI_COMM_WORLD, &world_size);


  PageRankType * thrd_sum = (PageRankType * )malloc( world_size * sizeof(PageRankType));
  uintV * split_ptr_begin = (uintV *) malloc ( world_size* sizeof(uintV));
  for(int i=0; i< world_size; i++){
    thrd_sum[i] = 0;
    split_ptr_begin[i] = 0;
  }
  

  // generating split_ptr_begin
  int j=0;
  for (int i=1; i< workers; i++){
    split_ptr_begin[i] = split_ptr_begin[i-1]+thrd_load[j]; 
    j++;
  }

  t1.start();
  //Calculate Page Rank----------------------------
    
    v_start = split_ptr_begin[world_rank];
    v_end = split_ptr_begin[world_rank] + thrd_load[world_rank];
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
    }
    for (uintV i=0; i< n; i++){
      pr_next[i] = 0.0;
    }
    }
    printf("%d, %d, %lf\n",id, edge_count, time_b1);


    // Each process calculates partial sum for its vertx range 
    for(uintV i=v_start; i<v_end; i++ ){
      local_sum += pr_curr[i];
    }

  MPI_Reduce(&local_sum, &tot_sum, 1, PAGERANK_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD );



  time_taken = t1.stop();
  
  free(split_ptr_begin);        //free the malloced stuff
  free(thrd_sum );

  if(world_rank == 0){
  std::cout << "Sum of page rank : " << tot_sum << "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  }
  delete[] pr_curr;
  delete[] pr_next;
  delete[] pr_holdr;
}
int pageRankParallel_1(Graph &g, int workers, int world_rank,  int max_iters, uintV * thrd_load, uintV *thrd_edge ) {
  uintV n = g.n_;

  PageRankType *pr_curr = new PageRankType[n];
  PageRankType *pr_next = new PageRankType[n];
  PageRankType *pr_holdr = new PageRankType[n];
  PageRankType *pr_holdr_temp = new PageRankType[n];
  
  timer t1;
  double time_taken = 0.0;
  PageRankType tot_sum = 0;
  PageRankType local_sum = 0;


    timer  t_b1;
    double time_b1 = 0.0;
    uintE out_degree;
    uintV v;
    uintE edge_count =0;
    uintV v_start, v_end;

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
    //pr_holdr[i] = 0.0;
  }

  //int world_rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size = workers;
  int id = world_rank;
  //MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  PageRankType * thrd_sum = (PageRankType * )malloc( world_size * sizeof(PageRankType));
  uintV * split_ptr_begin = (uintV *) malloc ( world_size* sizeof(uintV));
  for(int i=0; i< world_size; i++){
    thrd_sum[i] = 0;
    split_ptr_begin[i] = 0;
  }

  // generating split_ptr_begin
  int j=0;
  for (int i=1; i< workers; i++){
    split_ptr_begin[i] = split_ptr_begin[i-1]+thrd_load[j]; 
    j++;
  }

  t1.start();
  //Calculate Page Rank----------------------------
  


    
    
    v_start = split_ptr_begin[id];
    v_end = split_ptr_begin[id] + thrd_load[id]; 
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
    t_b1.start();         //// --- synchronization phase 1 start ---
    
    //MPI_Reduce( pr_next, pr_holdr, n, PAGERANK_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( pr_next, pr_holdr_temp, n, PAGERANK_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD );
    
    MPI_Scatterv( pr_holdr_temp, thrd_load, split_ptr_begin, PAGERANK_MPI_TYPE, pr_holdr , thrd_load[id], PAGERANK_MPI_TYPE, 0, MPI_COMM_WORLD);

    time_b1 += t_b1.stop();   //// --- synchronization phase 1 end -----
    uintV v;
    for ( v = v_start; v < v_end; v++) {
      pr_next[v] = PAGE_RANK(pr_holdr[v - split_ptr_begin[id] ]);           //this eliminates above for loop

      // reset pr_curr for the next iteration
      pr_curr[v] = pr_next[v];
    }
    for (v =0; v< n; v++){
      pr_next[v] = 0.0;
    }
    
    }
    time_taken = t1.stop();
    printf("%d, %d, %lf\n",id, edge_count, time_b1);

    // Each process calculates partial sum for its vertx range 
    for(uintV i=v_start; i<v_end; i++ ){
      local_sum += pr_curr[i];
    }
    


  MPI_Reduce(&local_sum, &tot_sum, 1, PAGERANK_MPI_TYPE, MPI_SUM, 1, MPI_COMM_WORLD );



  time_taken = t1.stop();
  
  if(world_rank == 1){
  std::cout << "Sum of page rank : " << tot_sum << "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  }
  
  free(split_ptr_begin);      //free malloc stuff
  free(thrd_sum );

  delete[] pr_curr;
  delete[] pr_next;
  delete[] pr_holdr;
  delete[] pr_holdr_temp ;
}




int main(int argc, char *argv[]) {

   cxxopts::Options options("page_rank_push", "Calculate page_rank using serial and parallel execution");
     options.add_options("", {
                                 {"nIterations", "Maximum number of iterations", cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
                                 {"strategy", "Strategy to be used", cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                                 {"inputFile", "Input graph file path", cxxopts::value<std::string>()->default_value("/scratch/input_graphs/roadNet-CA")},
                             });
 
     auto cl_options = options.parse(argc, argv);
     uint strategy = cl_options["strategy"].as<uint>();
     uint max_iterations = cl_options["nIterations"].as<uint>();
     std::string input_file_path = cl_options["inputFile"].as<std::string>();

  
  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Assuming at least 2 processes for this task
  if (world_size < 2) {
    fprintf(stderr, "World size must be greater than 1 for this program\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  if (world_rank == 1){
  #ifdef USE_INT
     std::printf("Using INT\n");
  #else
     std::printf("Using FLOAT\n");
  #endif
  std::cout << std::fixed;
  std::printf("World size : %d\n", world_size);
  std::printf("Communication strategy : %d\n", strategy);
  std::printf("Iterations : %d\n", max_iterations);
  std::cout << "rank, num_edges, communication_time\n";
  }
  
  Graph g;
  g.readGraphFromBinary<int>(input_file_path);
  
  uintV n = g.n_;
  //uintV * thrd_load = (uintV *) calloc( world_size , sizeof(int));
  //uintV * thrd_edge = (uintV *) calloc(world_size , sizeof(int));
  uintV thrd_load[ world_size ];
  uintV thrd_edge[world_size ];
  for (int i =0; i<world_size; i++){
    thrd_load[i] = 0;
    thrd_edge[i] = 0;
  }
  get_load_4_equal_edge_pr(std::ref(g), thrd_load, thrd_edge, world_size, n );
  
  if(strategy == 1){
    pageRankParallel_1(g, world_size, world_rank , max_iterations, thrd_load, thrd_edge);
  }else{
    pageRankParallel_2(g, world_size, world_rank ,max_iterations, thrd_load, thrd_edge);
  }

  //free(thrd_load);
  //free(thrd_edge);

  MPI_Finalize();  
  return 0;
}
