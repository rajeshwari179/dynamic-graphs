#include<iostream>
#include<stdio.h>
#include <stdlib.h>
#include <cstring>
#include <dirent.h> 
#include <sys/time.h>
#include <omp.h>
#include <chrono>
#include"graph_methods.h"
#include"core_maintenance.h"

using namespace std;
using namespace std::chrono;

#define CLOCK_PER_MS 1000

void FindThreads(vector<params> param);
void DeleteThreads();
void InsertThreads();
void TravelInsert(vector<pair<int,int>> rootEdges);
void InsertRemove(int k, int r);
void GetKSuperiorEdges(params param);
void TravelDelete(vector<pair<int,int>> rootEdges);
void  DeleteRemove(int k, int r);
Graph graph;
newGraph newgraph;
vector<vector<pair<int,int>>> superiorEdges;


//delete superior edges from new graph and insert them into original graph
void InsertEdgesIntoGraph()
{
	//int superiorEdgesSize = superiorEdges.size();
	for(int i = 0;i < superiorEdges.size();i ++){
	    //k-superior edges list
		vector<pair<int,int> > kedges = superiorEdges[i];
		//int ksize = kedges.size();
		for(int k = 0;k < kedges.size();k ++){
			int u = kedges[k].first;
			int v = kedges[k].second;
			//delete the edge from the insertion edge set
			if(!newgraph.deleteEdge(u,v)){
				cout<<"did not delete in edge set"<<u<<","<<v<<endl;
				return;
			}
			//add the edge to the original graph
			if(!graph.addEdge(u,v) || !graph.addEdge(v,u)){
				cout<<"did not insert into original"<<u<<","<<v<<endl;
				return;
			}
		}
	}
}

//delete superior edges from original graph
void DeleteEdgesFromGraph()
{
	for(int i = 0;i < superiorEdges.size();i ++){
		vector<pair<int,int> > kedges = superiorEdges[i];
	//	int ksize = kedges.size();
		for(int k = 0;k <  kedges.size();k ++){
			int u = kedges[k].first;
			int v = kedges[k].second;
			if(!newgraph.deleteEdge(u,v)){
				cout<<"did not delete in edge set"<<u<<","<<v<<endl;
				return;
			}
			if(!graph.deleteEdge(u,v) || !graph.deleteEdge(v,u)){
				cout<<"did not delete "<<u<<","<<v<<endl;
				return;
			}
		}
	}
}

//create threads to find superior edges
void FindThreads(vector<params> param)
{	
	int num_threads = param.size();
	//omp_set_num_threads(num_threads);
	//each thread is assigned a specific index, allowing them to work on different sets of superior edges concurrently
	#pragma omp parallel for
    for (int i = 0; i < num_threads; i++) {
        GetKSuperiorEdges(param[i]);
    }

}

//create threads to call RemoveEdge
void DeleteThreads()
{
	int num_threads = superiorEdges.size();
	//omp_set_num_threads(num_threads);
	#pragma omp parallel for
	for (int i = 0; i < num_threads; ++i) {
	// Each thread executes TravelDelete function with the kth superior edges list
		TravelDelete(superiorEdges[i]);
}
	
}

//create threads to call InsertEdge
void InsertThreads()
{
	int num_threads = superiorEdges.size();
	//omp_set_num_threads(num_threads);
	#pragma omp parallel for
    for (int i = 0; i < num_threads; ++i) {
    // Each thread executes TravelInsert function with the kth superior edges list
        TravelInsert(superiorEdges[i]);
    }

}

/*
given a core number k, find the k-superiorEdges in newgraph(the inserted/deleted edges) and add to superiorEdges
The index variable is used to determine which superior edge set to append the found superior edges to.
Each thread may be assigned a specific index, allowing them to work on different sets of superior edges concurrently.
*/
//param = k, i
void GetKSuperiorEdges(params param)
{
    //params par = *(params*)param;
	params par = param;
    int k = par.k;
    //index of k^th superior edge set in the list
    int index = par.i; 
    vector<int> rootvertex;
    vector<bool> IsPicked;
    //int verNumber = 0;
    //collate list of all vertices with core number k
    for(int i=0;i<newgraph.NumV;i++){
        if(newgraph.vertices[i].core == k){
            rootvertex.push_back(newgraph.vertices[i].id);
        }
        IsPicked.push_back(false);
    }
    //verNumber = rootvertex.size();
    for(int i = 0;i < rootvertex.size();i++){
        int idu = rootvertex[i];
        //index of the vertex u in the graph
        int index_u = newgraph.index_map[idu];
        //is picked = did we find a neighbour of u with core number >=k
		//every vertex can be picked only once to update the core number
		//this ensures it exists only in one superior edge set and we can parallize
        if(!IsPicked[index_u]){
            vector<int> neighbours = newgraph.vertices[index_u].neighbours;
            //int uSize = neighbours.size();
            for(int j=0;j<neighbours.size();j++){
                int index_v = neighbours[j];
                int idv = newgraph.vertices[index_v].id;
                int corev = newgraph.vertices[index_v].core;
                if(corev >= k && !IsPicked[index_v]){
                    IsPicked[index_u]=true;
                    //if corev == k, set v picked
                    if(corev == k){
                        IsPicked[index_v]=true;
                    }
                    superiorEdges[index].push_back(make_pair(idu,idv));
                    break;
                }
            }
        }
    }
}

/**
*after inserting a superior edge set, search vertices whose core numbers are same as the roots and
*find vertices that will change cores
*executed by a thread
*sameCoreEdges:the roots of new inserted edges, for each root conducts the InsertOneEdge algorithm
*+ve DFS is conducted on vertices in KPT(r) which has core number k to find the set of vertices whose core numbers potentially change**/
//Algorithm 2 - For each k, a child process is assigned to find the vertices whose core number changes are casused by the insertion of
//the k-superior edge set
void TravelInsert(vector<pair<int,int> > rootEdges)
{
    //root edges (k-superior edge set Ek)
    //vector<pair<int,int> > rootEdges=*(vector<pair<int,int> >*)sameCoreEdges;
    stack<int>s;
    int k = 0;
    int r = 0;
    //int rootsize = rootEdges.size();
    //compute SD(v) for each vertex v in exPT of Ek - Lemma 3
    graph.computeInsertSD(rootEdges);
    for(int i=0;i<rootEdges.size();i++){
        int u = rootEdges[i].first;
        int v = rootEdges[i].second;
        //int coreu = graph.vertices[u].core;
        //int corev = graph.vertices[v].core;
        r = u;
        k = graph.vertices[u].core;
        //find the vertex with smaller core number
        if( graph.vertices[u].core > graph.vertices[v].core ){
            r = v;
            k = graph.vertices[v].core;
        }
        if(!(graph.vertices[r].visited && !graph.vertices[r].removed)){
            graph.vertices[r].visited = true;
            //compute CSD (because only if CSD(w) > core(u) then w can have its core number increased)
            if(!graph.vertices[r].csd){
                graph.computeCSD(r);
            }
			//Lemma 4 - once a vertex has been updated, its core number will not increase anymore
			//cd = potential of a vertex to increase its core number. if cd(v)<=k, its core number cannot increase
            if(graph.vertices[r].cd >= 0)
                graph.vertices[r].cd = graph.vertices[r].csd;
            //if r has been visited before, its cd should be updated not initialized
            else
                graph.vertices[r].cd += graph.vertices[r].csd;
            s.push(r);
            int cdr = graph.vertices[r].cd;
            while(!s.empty()){
                int v=s.top();
                s.pop();
               // int cdv = graph.vertices[v].cd;
                //if cdv > k, then core number of v can increase
                if(graph.vertices[v].cd > k){
                   // int sizev = graph.edges[v].size();
                    for(int j=0;j<graph.edges[v].size();j++){
                        int w = graph.edges[v][j];
                        //int corew = graph.vertices[w].core;
                        /*if(!graph.vertices[w].sd){
                            graph.computeSD(w);
                        }*/
                        //if the core number of w is k and SD(w)>k, then vw is a k superior edge
                        //int sdw = graph.vertices[w].sd;
                        if(graph.vertices[w].core == k && graph.vertices[w].sd > k &&(!graph.vertices[w].visited)){
                            s.push(w);
                            graph.vertices[w].visited = true;
                            if(!graph.vertices[w].csd){
                                graph.computeCSD(w);
                            }
                            graph.vertices[w].cd += graph.vertices[w].csd;
                        }
                    }
                }
				//-ve DFS to remove v, update cd values of other vertices with core number k
				//once all vertices are traversed, vertices that are visited but not removed get core number++.
                else
                if(!graph.vertices[v].removed){
                    InsertRemove(k,v);
                }
            }
        }
    }

}
//missing: for each vertex v in G do:
//if !removed and visited then:
//add v to Vk and return Vk


//Algorithm 3
void  InsertRemove(int k, int r)
{
	stack<int> s;
    s.push(r);
    graph.vertices[r].removed = true;
    while(!s.empty()){
        int v = s.top();
        s.pop();
        for(int i=0;i<graph.edges[v].size();i++){
            int w = graph.edges[v][i];
            //int corew = graph.vertices[w].core;
            if(graph.vertices[w].core == k){
                graph.vertices[w].cd--;
                //int cdw = graph.vertices[w].cd;
                if(graph.vertices[w].cd == k && !graph.vertices[w].removed){
                    s.push(w);
                    graph.vertices[w].removed = true;
                }
            }
        }
    }
}

/**
*after deleting a superior edge set, search vertices whose core numbers are same as the roots and
*find vertices that will change cores
*sameCoreEdges:the roots of deleted edges**/
void TravelDelete(vector<pair<int,int> > rootEdges)
{
    //vector<pair<int,int> > rootEdges=*(vector<pair<int,int> >*)sameCoreEdges;
    int rootsize = rootEdges.size();
    for(int i = 0;i < rootsize;i ++){
        int u1 = rootEdges[i].first;
        int u2 = rootEdges[i].second;
        //int coreu1 = graph.vertices[u1].core;
        int coreu2 = graph.vertices[u2].core;
        int r = u1;
        int k = graph.vertices[u1].core;
        //update root as vertex with smaller core number
        if(graph.vertices[u1].core > graph.vertices[u2].core){
            r=u2;
            k=graph.vertices[u2].core;
        }
        //if core numbers are not equal
        if(graph.vertices[u1].core != graph.vertices[u2].core){
            
            if(!graph.vertices[r].visited ){
                graph.vertices[r].visited = true;
                //compute superior degree if not computed
                if(!graph.vertices[r].sd){
                    graph.computeSD(r);
                }
                graph.vertices[r].cd = graph.vertices[r].sd ;
            }
            if(!graph.vertices[r].removed){
                //int cdr = graph.vertices[r].cd;
                //if CD < k, remove the node
                if(graph.vertices[r].cd < k){
                    DeleteRemove(k,r);
                }
            }
        }
        else{
            //remove the node if the core number cd is less than k
            if(!graph.vertices[u1].visited){
                graph.vertices[u1].visited = true;
                if(!graph.vertices[u1].sd){
                    graph.computeSD(u1);
                }
                graph.vertices[u1].cd = graph.vertices[u1].sd;
            }
            if(!graph.vertices[u1].removed){
               //int cdu1 = graph.vertices[u1].cd;
               if(graph.vertices[u1].cd < k){
                    DeleteRemove(k,u1);
               }
           }
           if(!graph.vertices[u2].visited){
                graph.vertices[u2].visited = true;
                if(!graph.vertices[u2].sd){
                    graph.computeSD(u2);
                }
                graph.vertices[u2].cd = graph.vertices[u2].sd;
            }
            if(!graph.vertices[u2].removed){
               // int cdu2 = graph.vertices[u2].cd;
               if(graph.vertices[u2].cd < k ){
                    DeleteRemove(k,u2);
               }
           }
        }
    }

}

//remove the root node if its cd < k
//Perform DFS and remove neighbours if their cd < k
void  DeleteRemove(int k, int r)
{
    stack<int> s;
    s.push(r);
    graph.vertices[r].removed = true;
    while(!s.empty()){
        int v = s.top();
        s.pop();
        for(int i=0;i<graph.edges[v].size();i++){
            //neighbours of v
            int w = graph.edges[v][i];
            //int corew = graph.vertices[w].core;
            if(graph.vertices[w].core == k){
                if(!graph.vertices[w].visited){
                    if(!graph.vertices[w].sd){
                        graph.computeSD(w);
                    }
                    //int sdw = graph.vertices[w].sd;
                    graph.vertices[w].cd += graph.vertices[w].sd;
                    graph.vertices[w].visited = true;
                }
                //decrement cd of neighbour as current node is being deleted
                graph.vertices[w].cd--;
                //int cdw = graph.vertices[w].cd;
                //if w's cd < k, add it to the stack to remove.
                if(graph.vertices[w].cd < k && !graph.vertices[w].removed){
                    graph.vertices[w].removed = true;
					s.push(w);
                }
            }
        }
    }
}


int main(int argc,char* argv[])
{
    if(argc != 4){
		cout<<"-s/-p graph_filename edge_filename"<<endl;
		return 0;
	}
	string fname=argv[2];
	string edge_file = argv[3];
	string graphfile = fname + ".txt";
	string edgefile = edge_file + ".txt";
	string delfile = fname + "_deletion.txt";
	string insertfile = fname + "_insertion.txt";
	ifstream fingraph(graphfile.data());
	ifstream finedge(edgefile.data());
	vector<pair<int,int> > allNewEdges;
	vector<int> allcorenums;
	struct timeval t_start,t_end,f_start,f_end; 
	double dur=0.0;
	{//open files, read graph and edge file and compute core
		int a,b;
		if (!fingraph){
			cout <<  "Could not open "  << graphfile << endl;
			return 0;
		}
		while(fingraph>>a>>b){
		//cout <<  "opening "  << graphfile << endl;
		//cout<<a<<" "<<b<<endl;
			graph.addEdge(a,b);
			graph.addEdge(b,a);
		}
		//compute core numbers
		graph.ComputeCoreNumbers();
		allcorenums = graph.GetCoreNumbers();
		//get new graph and set cores
		if (!finedge){
			cout <<  "Could not open "  << edgefile  << endl;
			return 0;
		}
		while(finedge>>a>>b){
			allNewEdges.push_back(make_pair(a,b));
		}
		//new edges to be inserted
		newgraph.Map_index(allNewEdges);
		newgraph.SetCores(allcorenums);
	}

    //parallel algorithm
	if(strcmp(argv[1],"-p") == 0)
	{ //delete the edges first and then insert them back
		string output = fname + "_pll_outtime.txt";
		ofstream fout(output.data(),ios::app);
		fout<<"Original Graph\tEdge File\tdeletion_time\tfind_time\tdelete_iterations\tinsertion_time\tfind_time\tinsert_iterations\n";
		fout<<fname<<"\t"<<edge_file<<"\t";
		//number of edges in new graph
		int EdgeSet = newgraph.GetNumEdges();
		//cout<<EdgeSet<<endl;
		{//delete the edges
			gettimeofday(&t_start, NULL); 
			double startall_t = ((double)t_start.tv_sec)*1000+(double)t_start.tv_usec/1000; 
			double findTime = 0.0;
			int iteration = 0;
			while(EdgeSet){
				iteration++;
				vector<params> param;
				gettimeofday(&f_start, NULL); 
				double start_t =((double)f_start.tv_sec)*1000+(double)f_start.tv_usec/1000;
				//get new root core numbers
				param = newgraph.GetRootVertices();
				//number of threads = param.size
				superiorEdges.resize(param.size());
				//each thread runs GetKSuperiorEdges - so this is the time we are counting
				FindThreads(param);
				gettimeofday(&f_end, NULL); 
				double end_t =((double)f_end.tv_sec)*1000+(double)f_end.tv_usec/1000;
				dur =(end_t-start_t)/CLOCK_PER_MS;
				findTime += dur;
				//delete superior edges from new graph and original graph
				DeleteEdgesFromGraph();
				
				{//deleting threads and update cores for new graph
					DeleteThreads();
					//decrease core numbers for vertices that are visited and removed
					graph.delCores();
					//reset all vertices and get cores again
					graph.resetVertices();
					allcorenums = graph.GetCoreNumbers();
					//set cores for new graph
					newgraph.SetCores(allcorenums);
					//get number of edges in new graph
					EdgeSet = newgraph.GetNumEdges();
					superiorEdges.clear();
				}
			}
			gettimeofday(&t_end, NULL); 
			long endall_t = ((double)t_end.tv_sec)*1000+(double)t_end.tv_usec/1000; 
			dur = (endall_t-startall_t)/CLOCK_PER_MS;
			//write time taken for entire delete and GetKSuperiorEdges to the file
			fout<<dur<<"\t"<<findTime<<"\t"<<iteration<<"\t";	
			//write to new core file
			graph.WriteCores(delfile);
			
		}
		
		//Algorithm 1
		{//Compute Core Numbers in Original Graph
			graph.ComputeCoreNumbers();
			allcorenums = graph.GetCoreNumbers();
		//Construct the new vertex/edge set - (V',E')
			newgraph.clear();
			//initializes a graph data structure by mapping vertex IDs to indices,
			//creating vertex objects for new vertices 
			//updating the total vertex count and adding new edges to the graph
			newgraph.Map_index(allNewEdges);
			newgraph.SetCores(allcorenums);				
		}

		{//insertion
			gettimeofday(&t_start, NULL); 
			double startall = ((double)t_start.tv_sec)*1000+(double)t_start.tv_usec/1000;
			//E' = inserted edge set
			EdgeSet = newgraph.GetNumEdges();
			double findTime = 0.0;
			int iteration = 0;
			while(EdgeSet){
				iteration++;
				vector<params> param;
				gettimeofday(&f_start, NULL); 
				double start = ((double)t_start.tv_sec)*1000+(double)t_start.tv_usec/1000;
				//get root core numbers
				param = newgraph.GetRootVertices();
				//param has all vertices that have a neighbour with higher core number, i.e. a superior edge
				superiorEdges.resize(param.size());
				//each thread runs GetKSuperiorEdges to build Ek - this is the time we are counting in findTime
				FindThreads(param);
				gettimeofday(&f_end, NULL); 
				double end =((double)f_end.tv_sec)*1000+(double)f_end.tv_usec/1000;
				dur = (end-start)/CLOCK_PER_MS;
				findTime += dur;
				//delete superior edges from new graph and insert them into original graph
				InsertEdgesIntoGraph();
				{//Parallel K-Superior Insert
					InsertThreads();
				//increase core numbers for vertices that are visited but not removed -> parallel add?
					graph.insCores();
				//reset SD, CSD, CD for all vertices
					graph.resetVertices();
				//Recompute Core Numbers
					allcorenums = graph.GetCoreNumbers();
					newgraph.SetCores(allcorenums);
				//Recompute the inserted eedge set - this would have reduced now
					EdgeSet = newgraph.GetNumEdges();
					superiorEdges.clear();
			    }
			}
			gettimeofday(&t_end, NULL); 
			double endall = ((double)t_end.tv_sec)*1000+(double)t_end.tv_usec/1000; 
			dur = (endall-startall)/CLOCK_PER_MS;
			fout<<dur<<"\t"<<findTime<<"\t"<<iteration<<"\n";
			
			//write to new core file
			graph.WriteCores(insertfile);								
		}
		fout.close();
	}
	
	//serial
		else if(strcmp(argv[1],"-s") == 0){
		string output = fname + "_outtime.txt";
		ofstream fout(output.data(),ios::app);
		fout<<"Original Graph"<<"\t"<<"Edge File"<<"\t"<<"Deletion Time"<< "\t"<<"Insertion Time"<<"\n";
		fout<<fname<<"\t"<<edge_file<<"\t";
		{//delete the new edges from the original graph (if they exist)
			/*clock_t start,end;
			start = clock();*/
			//gettimeofday(&f_start, NULL); 
			//double start = ((double)t_start.tv_sec)*1000+(double)t_start.tv_usec/1000;
			auto start = high_resolution_clock::now();
			graph.Deletion(allNewEdges);
			/*end=clock();
			double dur =(double)(end-start)/CLOCK_PER_MS;*/
			auto end = high_resolution_clock::now();
			auto dura = duration_cast<microseconds>(end-start);
			//dur = (end-start)/CLOCK_PER_MS;
			cout<<"Deletion Time (microseconds): "<<dura.count()<<endl;
			fout<<dura.count()<<"\t";
			//write to core file
			graph.WriteCores(delfile);	
		}		
		{//insert new edges back			
			/*clock_t start,end;
			start = clock();*/
			gettimeofday(&f_start, NULL); 
			auto start = high_resolution_clock::now();
			graph.Insertion(allNewEdges);
			/*end=clock();
			double dur =(double)(end-start)/CLOCK_PER_MS;*/
			//double end =((double)f_end.tv_sec)*1000+(double)f_end.tv_usec/1000;
			auto end = high_resolution_clock::now();
			auto dura = duration_cast<microseconds>(end-start);
			//dur = (end-start)/CLOCK_PER_MS;
			cout<<"Insertion Time (microseconds): "<<dura.count()<<endl;
			fout<<dura.count()<<"\n";
			//write to new core file
			graph.WriteCores(insertfile);
		}
		fout.close();
		}
		
	{//close the input files
		fingraph.close();
		finedge.close();
		
	}	
	//cout<<" finished in "<<dur<<" sec"<<endl;
	return 0;
}
