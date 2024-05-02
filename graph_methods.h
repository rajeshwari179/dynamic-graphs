#include<vector>
#include<map>
#include<stack>
#include<time.h>
#include<algorithm>
#include<fstream>

using namespace std;

struct vertex{
    //core number
    int core;
    //superior degree
    int sd;
    //constraint superior degree
    int csd;
    //used to evaluate the potential of a vertex to increase its core number
    int cd;
    bool visited;
    bool removed;
    vertex(){
        visited=false;
        removed=false;
        core=0;
        sd=0;
        csd=0;
        cd=0;
    }
};

class Graph
{

//Number of vertices and edges
public: int NumV;
        int NumEdges;
        //core number of all vertices
		vector<int> allcorenums;

vector<vertex> vertices;
//adjacency list of edges
vector<vector<int>> edges;

Graph()
{
	NumV = 0;
	NumEdges = 0;
}

//initialize graph with v vertices and e edges
void init(int vNum, int eNum){
    NumV=vNum;
    NumEdges=eNum;
    for(int i=0;i<NumV;i++){
        vertex temp=vertex();
        vertices.push_back(temp);
        edges.push_back(vector<int>(0));
    }
}

void clear(){
    vertices.clear();
    edges.clear();
    NumV = 0;
	NumEdges = 0;
}

//add edge to graph
 bool addEdge(int from, int to)
{
    if(from == to)
        return false;
	while(from+1 > NumV) {
		NumV++;
		vertex temp=vertex();
        vertices.push_back(temp);
        edges.push_back(vector<int>(0));
    }
    //check if edge is already there
	for (int i = 0; i < edges[from].size(); i++) {
		if (edges[from][i] == to)
			return false;
	}
	//store the edge
	edges[from].push_back(to);
	NumEdges++;
	return true;
}

//delete edge from graph
 bool deleteEdge(int from, int to){
    //cleanliness check
    if(from > NumV){
        return false;
    }
    vector <int>::iterator Iter;
    for(Iter=edges[from].begin();Iter!=edges[from].end();Iter++){
      if (*Iter== to){
       edges[from].erase(Iter);
       Iter = edges[from].begin();       
       NumEdges--;
       return true;
      }
    }
    //if there is no such edge
    if(Iter == edges[from].end()){
        return false;
    }
    return false;
}

//total number of vertices in the graph
 int GetNumV()
{
	return NumV;
}

//total number of edges in the graph
 int GetNumEdges()
{
	return NumEdges;
}

//compute core number for the graph, static method
//core number = largest k st there exists a k-core containing u
//k - core = subgraph of vertices v where each deg(v)>=k 
void ComputeCoreNumbers(){
    int a,b;
    vector<int> bin,vpos,vsorted,deg;
    //deg = holds degree of i^th node
    //bin = count of vertices of each degree -> count array to sort vertices
    //vsorted = sorted array of vertices based on their degrees.
    //vpos = offset of vertex u in the sorted array
    a = b = 0;
    int n = NumV;
    int maxdegree = 0;
    for(int u = 0;u < n;u++){
        vpos.push_back(0);
        vsorted.push_back(0);
        int usize = edges[u].size();
        deg.push_back(usize);
        //find maxdegree for count sort
        if(usize > maxdegree){
            maxdegree = usize;
        }
    }
    for(int k = 0;k <= maxdegree;k ++){
        bin.push_back(0);
    }
    for(int u = 0;u < n;u ++){
        bin[deg[u]]++;
    }
    //bin is updated with cumulative count - bin[k] now holds the offset of the first vertex of degree k
    int start  = 0;
    for(int k = 0;k <= maxdegree;k++){
        int num = bin[k];
        bin[k] = start;
        start += num;
    }
    //reorder vertices based on degree
    //pos is used to keep track of where each vertex is placed after reordering.
    //bin is incremented because the next vertex of the same degree should be placed after the current vertex.
    for(int u = 0;u < n ;u ++){
        vpos[u] = bin[deg[u]];
        vsorted[vpos[u]] = u;
        bin[deg[u]]++;
    }
    
    //shift the offsets by 1 index to the right
    for(int d = maxdegree;d>0;d--){
        bin[d]=bin[d-1];
    }
    //If u is not already in the position pw in vert, it swaps the positions of u and w in vert,
//and also updates their positions in pos
    //Increments the count of vertices with degree du in bin (bin[du]++), 
    //as u has been moved to the next position.
//Decrements the degree of u, as it lost a neighbor with higher degree during 
//the swapping process.

    bin[0] = 0;
    int degu,posu,posw,w;
    for(int i = 0;i < n;i++){
        int v = vsorted[i];
        for(int j=0;j<edges[v].size();j++){
            //swap positions if neighbour with higher degree is found
            //update core number of neighbour if bigger one is found
            int u = edges[v][j];
            if(deg[u] > deg[v]){
                degu = deg[u];
                posu = vpos[u]; 
                //posw = bin[deg[u]] - starting offset 
                posw = bin[degu];
                //w = vsorted[bin[deg[u]]] - position of the first vertex with deg[u] in the sorted array
                w = vsorted[posw];
                if(u != w){
                    vpos[u] = posw;
                    //vsorted[pos[u]] = w
                    vsorted[posu] = w;
                    vpos[w] = posu;
                    //vsorted[pos[w]] = u
                    vsorted[posw] = u;
                }
                //bin[deg[u]] is incremented, then deg[u] is decremented -> this becomes core number
                bin[degu]++;
                deg[u]--;
            }
        }
    }
    //core_numbers
    for(int i=0;i < n;i++){
        vertices[i].core = deg[i];
    }
    //core numbers
	allcorenums = deg;
}

vector<int> GetCoreNumbers(){
	return allcorenums;
}

//serial methods

//delete edges serially
void Deletion(vector<pair<int,int> > allNewEdges){//delete edges	
	for(int k=0;k < allNewEdges.size();k ++){
		int a = allNewEdges[k].first;
		int b = allNewEdges[k].second;
		if(deleteEdge(a,b) && deleteEdge(b,a)){
		    //calculate superior degree SD for the vertices
			computeSD();
			RemoveEdge(a,b);
			resetVertices();
		}
	}		
}


//insert new edges serially
void Insertion(vector<pair<int,int> > allNewEdges){
	for(int k=0;k < allNewEdges.size();k ++){
		int a = allNewEdges[k].first;
		int b = allNewEdges[k].second;
		if(addEdge(a,b) && addEdge(b,a)){
		    //calculate SD
			computeSD();
			//calculate CSD
			computeCSD();
			InsertEdge(a,b);
			//core number of all vertices which are visited and not removed is incremented
			for(int i=0;i<NumV;i++){						
				if((vertices[i].visited)&&(!vertices[i].removed)){
					vertices[i].core++;
				}
			}
			resetVertices();
		}
	}
}

//deal with one edge insertion 
void  InsertEdge(int u1, int u2)
{
    //root = the vertex with smaller core number (should be vertex with core number k)
    int r=u1;
    int coreu1=vertices[u1].core;
    int coreu2=vertices[u2].core;
    if(coreu1 > coreu2){
        r=u2;
    }
    //positive DFS on KPT(r)
    stack<int> s;
    s.push(r);
    int K=vertices[r].core;
    vertices[r].visited=true;
    //initialize cd(v) as CSD(v).
    vertices[r].cd=vertices[r].csd;
    while(!s.empty()){
        int v=s.top();
        s.pop();
        int cdv = vertices[v].cd;
        //if cd(v)<=k, its core number cannot increase
        if(cdv > K){
            for(int j=0;j<edges[v].size();j++){
                int w=edges[v][j];
                int corew = vertices[w].core;
                //if the neighbour is counted in the constrained superior degree
                int sdw = vertices[w].sd;
                if(corew == K && sdw > K &&(!vertices[w].visited)){
                    s.push(w);
                    vertices[w].visited = true;
                    vertices[w].cd += vertices[w].csd;
                }
            }
        }
        else{
            if(!vertices[v].removed){
                InsertRemoveSerial(K,v);
            }
        }
    }
}

//Recursively remove vertices and decrement cd.
void  InsertRemoveSerial(int k, int v )
{
    vertices[v].removed=true;
    for(int i=0;i<edges[v].size();i++){
        int w=edges[v][i];
        //int corew = vertices[w].core;
        if(vertices[w].core == k){
            vertices[w].cd--;
            //int cdw = vertices[w].cd;
            //if w also has core number k, it has to be removed.
            if(vertices[w].cd == k && !vertices[w].removed){
                InsertRemoveSerial(k,w);
            }
		}
    }
}

//Single edge deletion
void RemoveEdge(int u1, int u2)
{
    int r = u1;
    int coreu1 = vertices[u1].core;
    int coreu2 = vertices[u2].core;
    int k = coreu1;
    if(coreu1 > coreu2){ 
		r=u2;
        k=coreu2;
    }
    if(coreu1 != coreu2){
        vertices[r].visited = true;
        vertices[r].cd = vertices[r].sd;
        int cdr = vertices[r].cd;
        if(cdr < k){
            DeleteRemove(k,r);
        }
    }
    else{
       vertices[u1].visited = true;
       vertices[u1].cd = vertices[u1].sd;
       int cdu1 = vertices[u1].cd;
       if(cdu1 < k){
            DeleteRemove(k,u1);
       }
       vertices[u2].visited = true;
       vertices[u2].cd = vertices[u2].sd;
       int cdu2 = vertices[u2].cd;
       if(!vertices[u2].removed && cdu2 < k ){
            DeleteRemove(k,u2);
       }
    }
}

void  DeleteRemove(int k, int v)
{
    vertices[v].removed = true;
    vertices[v].core--;
    for(int i=0;i< edges[v].size();i++){
        int w= edges[v][i];
        //int corew =   vertices[w].core;
        if(vertices[w].core == k){
            if(!  vertices[w].visited){
                  vertices[w].cd += vertices[w].sd;
                  vertices[w].visited = true;
            }
              vertices[w].cd--;
            //int cdw = vertices[w].cd;
            if(vertices[w].cd < k && ! vertices[w].removed){
                DeleteRemove(k,w);
            }
        }
    }
}


//Superior degree SD(u) = no of neighbours st v is a neighbour of u in G' and core(v)>=core(u).
//Only superior neighbours of a vertex can change the vertex's core number
void computeSD(){
    for(int v=0;v<NumV;v++){
        for(int j=0;j<edges[v].size();j++)
        {
            int w=edges[v][j];
            //int corev=vertices[v].core;
            //int corew=vertices[w].core;
            if(vertices[w].core >= vertices[v].core){
                vertices[v].sd++;
            }
        }
    }
}


//Constraint superior degree CSD(u) = core(w)>core(u) or core(w)==core(u) and SD(w)>core(u)
void computeCSD(){
    for(int v = 0;v < NumV;v++){
        for(int j=0;j<edges[v].size();j++)
        {
            int w=edges[v][j];
            //int corev=vertices[v].core;
            //int corew=vertices[w].core;
            //int sdw=vertices[w].sd;
            if(vertices[w].core > vertices[v].core ||
               (vertices[w].core == vertices[v].core && vertices[w].sd> vertices[v].core)){
                    vertices[v].csd++;
            }
        }
    }
}


//parallel methods


//decrease core numbers for vertices that are visited and removed
void delCores(){
	for(int i=0;i<NumV;i++){
		if((vertices[i].visited)&&(vertices[i].removed)){
			vertices[i].core--;
			allcorenums[i]--;
		}
	}
}

//increase core numbers for vertices that are visited but not removed
void insCores(){
	for(int i=0;i<NumV;i++){
		if((vertices[i].visited)&&(!vertices[i].removed)){
			vertices[i].core++;
			allcorenums[i]--;
		}
	}
}
//vectorize

/*void insCores(){
    std::vector<int> coresToAdd(NumV, 0);
	for(int i=0;i<NumV;i++){
		if((vertices[i].visited)&&(!vertices[i].removed)){
            coresToAdd[i] = 1;
		}
	}
    //vectorized addition
    #pragma omp simd
    for(int i=0;i<NumV;i++){
        vertices[i].core += coresToAdd[i];
        allcorenums[i] -= coresToAdd[i];
    }
}*/



//given a superior edge set, compute the (superior degree) for vertices in the EXPTree
//Used in Algorithm 2
//Build K-path tree rooted at r
 void computeInsertSD(vector<pair<int,int>> superioredges){
    map<int,bool> visited;	
    for(int i = 0;i < superioredges.size();i++){
        stack<int> s;
        int a = superioredges[i].first;
        int b = superioredges[i].second;
       // int sizea = edges[a].size();
        //int sizeb = edges[b].size();
        //int corea = vertices[a].core;
        //int coreb = vertices[b].core;
        //root vertex is the one with smaller core number - because this would have more superior edges
        int r = a;
        int core_r = vertices[a].core;
        int edges_r = edges[a].size();
        if(vertices[a].core > vertices[b].core){
            r = b;
            core_r = vertices[b].core;
            edges_r = edges[b].size();
        }
        //visit 
        //if(!visited[r] && !removed[r])
		if(!visited[r]){
			s.push(r);
			visited[r] = true;
			//check how many neighbours of r have higher core number and update SD of r
			for(int j=0;j<edges_r;j++)
			{
				int w = edges[r][j];
				int core_w = vertices[w].core;
				if(core_w >= core_r){
					vertices[r].sd++;
				}
			}
		   while(!s.empty()){
				int v=s.top();
				s.pop();
				//int sizev = edges[v].size();
				//int corev = vertices[v].core;
				//v's neighbours - check how many have same core number and visit them
				for(int j=0;j<edges[v].size();j++){
					int p = edges[v][j];
					if(vertices[p].core == vertices[v].core && !visited[p]){
						s.push(p);
						//vertices[p].visited = true;
						visited[p] = true;
						//compute superior degree if it hasn't been computed
						if(!vertices[p].sd){
							//int sizep = edges[p].size();
							//int corep = vertices[p].core;
							for(int k = 0;k < edges[p].size();k ++){
								int w = edges[p][k];
								//int corew = vertices[w].core;
								if(vertices[w].core >= vertices[p].core){
									vertices[p].sd++;
								}
							}
						}
					}
				}
			}
		}
	}

}


 void resetVertices(){
    for(int v = 0;v < NumV;v++){
        vertices[v].visited=false;
        vertices[v].removed=false;
        vertices[v].sd=0;
        vertices[v].csd=0;
        vertices[v].cd=0;
    }
}

//compute SD (Superior Degree) for a vertex v
 void computeSD(int v){
    //int corev = vertices[v].core;
    //int sizev = edges[v].size();
    for(int j=0;j<edges[v].size();j++)
    {
        int w = edges[v][j];
        //int corew = vertices[w].core;
        if(vertices[w].core >= vertices[v].core){
            vertices[v].sd++;
        }
    }
}


//compute CSD (Constraint Superior Degree) for a vertex v
 void computeCSD(int v){
   // int corev = vertices[v].core;
    //int sizev = edges[v].size();
    for(int j=0;j< edges[v].size();j++){
        int w = edges[v][j];
        //int corew = vertices[w].core;
        //int sdw = vertices[w].sd;
        if(vertices[w].core > vertices[v].core || (vertices[w].core == vertices[v].core && vertices[w].sd > vertices[v].core)){
            vertices[v].csd++;
        }
    }
}

//write core numbers to file
void WriteCores(string corefile){
	ofstream fcore(corefile.data());
	vector<pair<int,int> > vertex_cores;
	for(int i=0;i<NumV;i++){
		vertex_cores.push_back(make_pair(i,vertices[i].core));
	}
	
	auto comp = [](const pair<int, int>& a, const pair<int, int>& b) {
        return a.first < b.first;
    };
	
	sort(vertex_cores.begin(),vertex_cores.end(),comp);
	
	//write vertex ID, core number
	for(int i=0;i<vertex_cores.size();i++){
		fcore<<vertex_cores[i].first<<","<<vertex_cores[i].second<<endl;
	}
	fcore.close();
}




};

