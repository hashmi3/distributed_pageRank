# Distributed PageRank using OpenMPI

This repo consists of distributed programs built using cumulative message passing feature offered by OpenMPI library.

## Page Rank using OpenMPI
Implemented [Page rank](https://en.wikipedia.org/wiki/PageRank) algorithm over a compute cluster in a distributed manner. The work is distributed among `P` processes. For simplicity, every process will read the entire input graph and create its own graph data structure. For questions `3` and `4`, we will use **edge decomposition** strategy to distribute the work among the processes.


The PageRank pseudocode is:

    for each process P in parallel {
        communication_time = 0.0;
        for(i=0; i<max_iterations; i++) {
            for each vertex 'u' allocated to P {  // Loop 1
                edges_processed += outDegree(u)  // used in output validation
                for vertex 'v' in outNeighbor(u){
                    next_page_rank[v] += (current_page_rank[u]/outdegree[u])
                }
            }

            //Only the time spent on synchronization phase 1 is used for calculating communication time.
            // --- synchronization phase 1 start ---
            timer1.start();
            for each vertex 'u' allocated to P, aggregate (i.e., sum up) the value of next_page_rank[u] from all processes
            communication_time += timer1.stop();
            // --- synchronization phase 1 end -----
    
            for each vertex 'v' allocated to P {  // Loop 2
                new_value = PAGERANK(next_page_rank[v])
                current_page_rank[v] = new_value
            }
            Reset next_page_rank[v] to 0 for all vertices
        }
    }
    local_sum = 0
    for each vertex 'v' allocated to P {  // Loop 3
        local_sum += current_page_rank[v]
    }
    // --- synchronization phase 2 start ---
    if(P is root process){
        global_sum = Aggregated value of local_sum from all the processes using MPI_Reduce
        // print process statistics and other results
    }
    else{
        // print process statistics.
    }
    // --- synchronization phase 2 end ---
    
    
Two different strategies were adapted for calculating page rank. 
*   Strategy 1 -- `next_page_rank` is calculated at root node and later on dispersed to appropriate node.
*   Strategy 2 -- Here appropriate values of `next_page_rank` according to vertex distribution is calculated by each node 


### 1\. PageRank - Reduce and Scatter \[25 points\]

 Strategy 1 take advantage of [MPI\_Reduce()](https://www.open-mpi.org/doc/current/man3/MPI_Reduce.3.php) and [MPI\_Scatterv()](https://www.open-mpi.org/doc/current/man3/MPI_Scatterv.3.php) for communication among distributed processes. In this strategy,

*   The root process sums up the `next_page_rank` value for every vertex in the graph using [MPI\_Reduce()](https://www.open-mpi.org/doc/current/man3/MPI_Reduce.3.php).
*   Then, the root process sends the aggregated `next_page_rank` value of each vertex to its appropriate process. For example, if vertex `v` is allocated to process `P1`, the aggregated `next_page_rank[v]` is sent only to `P1`.

Implemented as cli argument `--strategy 1`.

### 2\. PageRank - Direct Reduce \[25 points\]

In previous strategy, the `next_page_rank` value was summed up for each vertex on the root process and it dispersed this aggregated value to the appropriate process. Instead, in this strategy,

*   Each process sums up the `next_page_rank` value, for every vertex allocated to that process, using [MPI\_Reduce()](https://www.open-mpi.org/doc/current/man3/MPI_Reduce.3.php).

Implement as cli argument `--strategy 2`.
