#include "core/graph.h"
#include "core/utils.h"
#include <future>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <stdlib.h>
#include <stdio.h>
#include "utilit.h"
#include <mpi.h>

void calc_tri_count_thrd_2(Graph &g, int id, uintV v_start, uintV v_end,  uint  num_edges, long &root_count , int world_size){
  timer t1;
  double time_taken = 0.0;
  long count = 0;
  long * count_list;
  uintE out_edge_count = 0;
  uint k=0;
  
  t1.start();
  // Process each edge <u,v>
  for (uintV u = v_start; u < v_end; u++) {
    // For each outNeighbor v, find the intersection of inNeighbor(u) and
    // outNeighbor(v)
    uintE out_degree = g.vertices_[u].getOutDegree();
    out_edge_count += out_degree;
    for (uintE i = 0; i < out_degree; i++) {
      uintV v = g.vertices_[u].getOutNeighbor(i);
      count += countTriangles(g.vertices_[u].getInNeighbors(),
                                       g.vertices_[u].getInDegree(),
                                       g.vertices_[v].getOutNeighbors(),
                                       g.vertices_[v].getOutDegree(), u, v);
    }
  }
  
  time_taken = t1.stop();

  if (id == 0){
    count_list = (long*) malloc( world_size* sizeof(long));

  }
  
  MPI_Gather(&count, 1, MPI_LONG, count_list, 1, MPI_LONG, 0, MPI_COMM_WORLD );

  if(id ==0){
    for(k =0; k< world_size; k++){
      root_count += count_list[k];
    }
  }

  printf("%d, %u, %ld, %lf\n",id,  num_edges, count, time_taken);
  if(id == 0){
    free(count_list );
  }
} 

void calc_tri_2(Graph &g, int id, uintV v_start, uintV v_end,  uint  num_edges, long &root_count , int world_size){
  
  timer t1;
  double time_taken = 0.0;
  long count = 0;
  long * count_list;
  long holdr;
  uintE out_edge_count = 0;
  uint k=0;
  
  t1.start();
  // Process each edge <u,v>
  for (uintV u = v_start; u < v_end; u++) {
    // For each outNeighbor v, find the intersection of inNeighbor(u) and
    // outNeighbor(v)
    uintE out_degree = g.vertices_[u].getOutDegree();
    out_edge_count += out_degree;
    for (uintE i = 0; i < out_degree; i++) {
      uintV v = g.vertices_[u].getOutNeighbor(i);
      count += countTriangles(g.vertices_[u].getInNeighbors(),
                                       g.vertices_[u].getInDegree(),
                                       g.vertices_[v].getOutNeighbors(),
                                       g.vertices_[v].getOutDegree(), u, v);
    }
  }
  
  time_taken = t1.stop();

  
  MPI_Reduce(&count, &root_count, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
  //MPI_Reduce(&count, &holdr, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
  //root_count =  holdr;
  if(id ==0)
   printf(">>> Total triangles : %ld\n", root_count);
  printf("%d, %u, %ld, %lf\n",id,  num_edges, count, time_taken);

}


int triangleCountParallel_st_2(Graph &g, int n_workers) {
  
  uintV n = g.n_;
  double time_taken = 0.0;
  timer t1, partition1;
  double partition_time = 0.0;

  int * thrd_load = (int *) calloc(n_workers, sizeof(int));
  int * split_ptr_begin = (int*) calloc(n_workers, sizeof(int));  //list split begin
  int * thrd_edge = (int *) calloc(n_workers, sizeof(int));       //list split end
  
  //int thrd_load[n_workers];
  //int split_ptr_begin[n_workers];
  //int thrd_edge[n_workers];

  long triangle_count = 0;            
  
  long root_count = 0;
  uintE tot_edg = 0;
  uintE temp = 0;
  
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
  
  partition1.start();
  get_load_4_equal_edge(std::ref(g),thrd_load, thrd_edge, n_workers, n);
  partition_time = partition1.stop();

  int j=0;
  for (int i=1; i< n_workers; i++){
    split_ptr_begin[i] = split_ptr_begin[i-1]+thrd_load[j]; 
    j++;
  }
  uintV cur = 0;
  if (world_rank == 0)
    printf("rank, edges, triangle_count, communication_time\n");
  // The outNghs and inNghs for a given vertex are already sorted

  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  
  
  t1.start();
  

    calc_tri_2( std::ref(g), world_rank, split_ptr_begin[world_rank], split_ptr_begin[world_rank] +thrd_load[world_rank], thrd_edge[world_rank], triangle_count, world_size  );
      

  time_taken = t1.stop();
  

  free(thrd_load);
  free(thrd_edge);
  free(split_ptr_begin );

  if(world_rank != 0)
    return 0;    
  // Print the overall statistics
  std::cout << "Number of triangles : " << triangle_count << "\n";
  std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
  //printf("Partitioning time (in seconds) : %lf\n", partition_time);
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << time_taken << "\n";
  return 0;
  


}

int triangleCountParallel_st_1(Graph &g, int n_workers) {
  uintV n = g.n_;
  double time_taken = 0.0;
  timer t1, partition1;
  double partition_time = 0.0;

  int * thrd_load = (int *) calloc(n_workers, sizeof(int));
  int * split_ptr_begin = (int*) calloc(n_workers, sizeof(int));  //list split begin
  int * thrd_edge = (int *) calloc(n_workers, sizeof(int));       //list split end
  
  //int thrd_load[n_workers];
  //int split_ptr_begin[n_workers];
  //int thrd_edge[n_workers];

  long triangle_count = 0;            
  
  long root_count = 0;
  uintE tot_edg = 0;
  uintE temp = 0;
  
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
  
  partition1.start();
  get_load_4_equal_edge(std::ref(g),thrd_load, thrd_edge, n_workers, n);
  partition_time = partition1.stop();

  int j=0;
  for (int i=1; i< n_workers; i++){
    split_ptr_begin[i] = split_ptr_begin[i-1]+thrd_load[j]; 
    j++;
  }
  uintV cur = 0;
  if (world_rank == 0)
    printf("rank, edges, triangle_count, communication_time\n");
  // The outNghs and inNghs for a given vertex are already sorted

  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  
  
  t1.start();
  

    calc_tri_count_thrd_2( std::ref(g), world_rank, split_ptr_begin[world_rank], split_ptr_begin[world_rank] +thrd_load[world_rank], thrd_edge[world_rank], triangle_count, world_size  );
      

  time_taken = t1.stop();
  
  

  free(thrd_load);
  free(thrd_edge);
  free(split_ptr_begin );

  if(world_rank != 0){
    return 0;
  }    
  
  // Print the overall statistics
  std::cout << "Number of triangles : " << triangle_count << "\n";
  std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
  //printf("Partitioning time (in seconds) : %lf\n", partition_time);
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << time_taken << "\n";
  return 0;
}

int main(int argc, char *argv[]) {
  cxxopts::Options options("triangle_counting_serial", "Count the number of triangles using serial and  parallel execution");
     options.add_options("custom", {
                                       {"strategy", "Strategy to be used", cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                                       {"inputFile", "Input graph file path", cxxopts::value<std::string>()->default_value("/scratch/input_graphs/roadNet-CA")},
                                   });
 
     auto cl_options = options.parse(argc, argv);
     uint strategy = cl_options["strategy"].as<uint>();
     std::string input_file_path = cl_options["inputFile"].as<std::string>();
 
  std::cout << std::fixed;
  
  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Assuming at least 2 processes for this task
  /*if (world_size < 2) {
    fprintf(stderr, "World size must be greater than 1 for this program\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }*/


  if (world_rank == 1){
  std::cout << "World size : " << world_size << "\n";
  std::cout << "Communication strategy : "<< strategy << "\n";
  }

  Graph g;
  g.readGraphFromBinary<int>(input_file_path);
  

  //triangleCountParallel_2(g, 4);      //Testing
  //triangleCountParallel_2(g, n_workers );    

  if (strategy == 1){
    triangleCountParallel_st_1(g, world_size);
    //triangleCountParallel_st_1(g, 4);
  } else{
    triangleCountParallel_st_2(g, world_size );
    //triangleCountParallel_st_2(g, n_workers);
  } 

  //printf("Done with task ! from process: %d \n", world_rank);
  // Finalize the MPI environment. No more MPI calls can be made after this
  
  MPI_Finalize();  

  return 0;
}
