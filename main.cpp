#include<iostream>
#include<stdio.h>
#include <stdlib.h>
#include <cstring>
#include <dirent.h> 
#include <sys/time.h>
#include"graph_methods.h"
#include"core_maintenance.h"

using namespace std;

#define CLOCK_PER_MS 1000

void EdgeThreads(vector<params> param);
void DeleteThread();
void InsertThread();
void GetKSuperiorEdges(void* param);
Graph graph;
newGraph newgraph;
vector<vector<pair<int,int>>> superiorEdges;


//delete superior edges from new graph and insert them into original graph
void InsertEdgesIntoGraph()
{
	int superiorEdgesSize = superiorEdges.size();
	for(int i = 0;i < superiorEdgesSize;i ++){
		vector<pair<int,int> > kedges = superiorEdges[i];
		int ksize = kedges.size();
		for(int k = 0;k < ksize;k ++){
			int u = kedges[k].first;
			int v = kedges[k].second;
			if(!newgraph.deleteEdge(u,v)){
				cout<<"did not delete "<<u<<","<<v<<endl;
				return;
			}
			if(!graph.addEdge(u,v) || !graph.addEdge(v,u)){
				cout<<"did not insert "<<u<<","<<v<<endl;
				return;
			}
		}
	}
}

//delete superior edges from new graph and original graph
void DeleteEdgesFromGraph()
{
	for(int i = 0;i < superiorEdges.size();i ++){
		vector<pair<int,int> > kedges = superiorEdges[i];
		int ksize = kedges.size();
		for(int k = 0;k < ksize;k ++){
			int u = kedges[k].first;
			int v = kedges[k].second;
			if(!graph.deleteEdge(u,v) || !graph.deleteEdge(v,u)){
				cout<<"did not delete "<<u<<","<<v<<endl;
				return;
			}
		}
	}
}

//PENDING: create threads to find superior edges
void EdgeThreads(vector<params> param)
{

}

//PENDING: create threads to call RemoveEdge
void DeleteThreads()
{
	
}

//PENDING: create threads to call InsertEdge
void InsertThreads()
{

}

/*
given a core number k, find the k-superiorEdges in newgraph(the inserted/deleted edges) and add to superiorEdges
The index variable is used to determine which superior edge set to append the found superior edges to.
Each thread may be assigned a specific index, allowing them to work on different sets of superior edges concurrently.
*/
//param = k, i
void GetKSuperiorEdges(void* param)
{
    params par = *(params*)param;
    int k = par.k;
    int index = par.i; 
    vector<int> rootvertex;
    vector<bool> IsPicked;
    int verNumber = 0;
    //collate list of all vertices with core number k
    for(int i=0;i<newgraph.NumV;i++){
        if(newgraph.vertices[i].core == k){
            rootvertex.push_back(newgraph.vertices[i].id);
        }
        IsPicked.push_back(false);
    }
    verNumber = rootvertex.size();
    for(int i = 0;i < verNumber;i++){
        int idu = rootvertex[i];
        //index of the vertex u in the graph
        int index_u = newgraph.index_map[idu];
        //is picked = did we find a neighbour of u with core number >=k
        if(!IsPicked[index_u]){
            vector<int> adjacent = newgraph.vertices[index_u].adjacent;
            int uSize = adjacent.size();
            for(int j=0;j<uSize;j++){
                int index_v = adjacent[j];
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



int main(int argc,char* argv[])
{
    if(argc != 4){
		cout<<"-s graph_filename edge_filename"<<endl;
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
	double dur;
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
	
	//serial
		string output = fname + "_outtime.txt";
		ofstream fout(output.data(),ios::app);
		fout<<fname<<"\t"<<edge_file<<"\t";
		{//delete the new edges from the original graph (if they exist)
			clock_t start,end;
			start = clock();
			graph.Deletion(allNewEdges);
			end=clock();
			double dur =(double)(end-start)/CLOCK_PER_MS;
			cout<<"Deletion Time (ms): "<<(double)end-start<<endl;
			fout<<dur<<"\t";
			//write to core file
			graph.WriteCores(delfile);	
		}		
		{//insert new edges back			
			clock_t start,end;
			start = clock();
			graph.Insertion(allNewEdges);
			end=clock();
			double dur =(double)(end-start)/CLOCK_PER_MS;
			cout<<"Insertion Time (ms): "<<dur<<endl;
			fout<<dur<<"\n";
			//write to new core file
			graph.WriteCores(insertfile);
		}
		fout.close();
	
		
	{//close the input files
		fingraph.close();
		finedge.close();
		
	}	
	//cout<<" finished in "<<dur<<" sec"<<endl;
	return 0;
}
