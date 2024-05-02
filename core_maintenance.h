#include<vector>
#include<map>
#include<iostream>
#include<stack>
#include<algorithm>
#include<time.h>
using namespace std;


struct newvertex{
    int id;
    int core;

    //adjacency list
    vector<int> neighbours;
    newvertex(){
        id = core = 0;
        neighbours.resize(0);
    }

};

struct params{
    int k;
    int i;
    params(){
        k = i = 0;
    }
    params(int K,int index){
        k = K;
        i = index;
    }
};

class newGraph
{
public: int NumV;
        int NumEdges;
        map<int,int> index_map;
        vector<newvertex> vertices;
		
public: newGraph()
{
	NumV = NumEdges = 0;
}

public: void init(int vNum, int eNum){
    NumV=vNum;
    NumEdges=eNum;
    for(int i=0;i<NumV;i++){
        newvertex temp=newvertex();
        vertices.push_back(temp);
	}
}

public: void clear(){
    vertices.clear();
    index_map.clear();
}

 //initializes a graph data structure by mapping vertex IDs to indices,
 //creating vertex objects for new vertices
 //updating the total vertex count, and adding edges to the graph
public: void Map_index(vector<pair<int,int> > allNewEdges){
   // int from, to, index;
int from = 0;
  int  to = 0;
   int index = 0;
    for(size_t i=0;i<allNewEdges.size();i++){
         from = allNewEdges[i].first;
         to = allNewEdges[i].second;
        //if either vertex is not in the map, add and create new vertex object
        if(index_map.find(from) == index_map.end()){
            index_map[from] = index;
            newvertex temp=newvertex();
            temp.id = from;
            vertices.push_back(temp);
            index++;
        }
        if(index_map.find(to) == index_map.end()){
            index_map[to] = index;
            newvertex temp=newvertex();
            temp.id = to;
            vertices.push_back(temp);
            index++;
        }
    }
    //total number of vertices
    NumV = index_map.size();
    
    //add edges to graph
	for(int i=0;i<allNewEdges.size();i++){
		//int a = allNewEdges[i].first;
		//int b = allNewEdges[i].second;
		addEdge(allNewEdges[i].first,allNewEdges[i].second);
		addEdge(allNewEdges[i].second,allNewEdges[i].first);
	}
}

public: bool addEdge(int from, int to)
{
	if(from == to )
        return false;
    //find the index of vertex from
	int index_f = index_map[from];
	int index_t = index_map[to];
	
	int len = vertices[index_f].neighbours.size();
	for (int i = 0; i < len; i++) {
		if (vertices[index_f].neighbours[i] == index_t)
			return false;
	}

	//store the edge
	vertices[index_f].neighbours.push_back(index_t);
	NumEdges++;
	return true;
}

public: bool deleteEdge(int from, int to){
    bool flag = true;
    int index_f = index_map[from];
	int index_t = index_map[to];
    if(index_f > NumV){
        flag = false;
    }
    //delete from vertex
    vector <int>::iterator Iter = std::find(vertices[index_f].neighbours.begin(),vertices[index_f].neighbours.end(), index_t);
    if(Iter != vertices[index_f].neighbours.end()){
        vertices[index_f].neighbours.erase(Iter);
        NumEdges--;
    }
    else{
        flag = false;
    }
    //delete to vertex
    vector<int>::iterator it1 = std::find(vertices[index_t].neighbours.begin(),vertices[index_t].neighbours.end(), index_f);
    if(it1 != vertices[index_t].neighbours.end()){
        vertices[index_t].neighbours.erase(it1);
        NumEdges--;
        flag = true;
    }
    else{
        flag = false;
    }
    return flag;
}

public: int GetNumV()
{
	return NumV;
}

public: int GetNumEdges()
{
	return NumEdges;
}

public: int GetVertexId(int u)
{
	return vertices[u].id;
}

int GetVertexCore(int u)
{
	return vertices[u].core;
}

vector<int> GetVertexAdj(int u)
{
	return vertices[u].neighbours;
}

//set cores for new graph, if its a new vertex, core is 0
public: void SetCores(vector<int> allcorenumbers){
	for(int u=0;u<NumV;u++){
		int idu = vertices[u].id;
		if(idu+1 > allcorenumbers.size()){
		//if(idu >= allcorenumbers.size()){
			vertices[u].core = 0;
		}else{
			vertices[u].core = allcorenumbers[idu];
		}
	}
}

//Used in Algorithm 1
//get root core numbers and first vertex with that core number
//find vertices with superior edges
public: vector<params> GetRootVertices(){
	vector<int> cores;
	vector<params> parameters;
	int vertex_index = 0;
	for(int u=0;u<NumV;u++){
	//	int coreu = vertices[u].core;
		int root_core = vertices[u].core;
		//vertex u has a neighbor with larger core number
		bool flag = false;
		for(int i=0;i<vertices[u].neighbours.size();i++){
			int v = vertices[u].neighbours[i];
			//int corev = vertices[v].core;
			if(vertices[v].core >= vertices[u].core){
				flag = true;
				break;
			}
		}
        //if u connects to a superior edge in E' and core(u) is not in C
		//if u has a neighbour with higher core number, add u's core number to root_cores list if its not there
		//(basically find vertices with superior edges)
		if(flag){
			vector<int>::iterator iter=std::find(cores.begin(),cores.end(),root_core);
			// if the core is not in root_cores, add it
			if(iter == cores.end()){
				cores.push_back(root_core);
				parameters.push_back(params(root_core,vertex_index));
				vertex_index++;
			}
		}
	}
	return parameters;
}
};
