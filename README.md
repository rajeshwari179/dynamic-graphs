# parallel-algorithms-for-dynamic-graphs
CSE6230-High Performance Parallel Computing Project

### Overview

This repository contains a program that implements methods for computing and maintaining core numbers of vertices in a graph, before/after deleting and inserting edges.

### To Build

```bash
make
```

### To Run

```bash
./mykcores -s <original graph filename> <list of edges to be inserted>
```
### Data Format

Both files should be in .txt format and the data should be in the following format:
<from vertex 1> <to vertex 1>
<from vertex 2> <to vertex 2>
.
.
And so on.

### Output Files
Three output files are generated:
```bash
1.	<original graph filename>_deletion.txt: This file contains the list of each vertex and its core number after the edges in <list of edges to be inserted> have been deleted.
2.	<original graph filename>_insertion.txt: This file contains the list of each vertex and its core number after the edges in <list of edges to be inserted> have been inserted.
3.	<original graph filename>_outtime.txt: This file contains the time (in ms) taken to compute the core numbers in the following format:
<original graph filename> <list of edges to be inserted> <time taken to delete edges and recompute core numbers> <time taken to insert edges and recompute core numbers>
```

### Code Files
graph_methods.h – this header file contains methods to initialize the graph, add/remove edges, compute the core number of vertices, compute the superior degree and constraint superior degree of vertices and write the core numbers to an output file.

core_maintenance.h – this header file contains methods to initialize the same graph with a structure that can be used to find vertices with superior edges efficiently after new edges have been inserted.

main.cpp – has methods to find the k-superior edges and insert/delete superior edges into the graph.
