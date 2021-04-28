#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

#include<iostream>
#include<string>
#include<vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>
#include <stack>
#include <limits>
#include <utility> 
#include <algorithm>
#include <string>
// include more libraries here if you need to

using namespace std; // the standard namespace are here just in case.

/*
	the vertex class
*/
template <typename T>
class vertex {

public:
	int id;
	T weight;

	vertex(int x, T y) : id(x), weight(y) {}

	// add more functions here if you need to
};

/*
	the graph class
*/
template <typename T>
class directed_graph {

private:
	vector<T> vertices_weights; //Stores the weight of the matrix per id using index
	vector<vector<T>> adj_matrix; //Stores the edges weight using matrix coordinates [i],[j] = weight where i and j is a vertex id 
	//You will need to add some data members here
	//to actually represent the graph internally,
	//and keep track of whatever you need to.

public:

	directed_graph(); //A constructor for directed_graph. The graph should start empty.
	~directed_graph(); //A destructor. Depending on how you do things, this may not be necessary.

	void increaseCapacity();

	bool contains(const int&) const; //Returns true if the graph contains the given vertex_id, false otherwise.
	bool adjacent(const int&, const int&) const; //Returns true if the first vertex is adjacent to the second, false otherwise.

	void add_vertex(const vertex<T>&); //Adds the passed in vertex to the graph (with no edges).
	void add_edge(const int&, const int&, const T&); //Adds a weighted edge from the first vertex to the second.

	void remove_vertex(const int&); //Removes the given vertex. Should also clear any incident edges.
	void remove_edge(const int&, const int&); //Removes the edge between the two vertices, if it exists.

	size_t in_degree(const int&) const; //Returns number of edges coming in to a vertex.
	size_t out_degree(const int&) const; //Returns the number of edges leaving a vertex.
	size_t degree(const int&) const; //Returns the degree of the vertex (both in edges and out edges).

	size_t num_vertices() const; //Returns the total number of vertices in the graph.
	size_t num_edges() const; //Returns the total number of edges in the graph.

	vector<vertex<T>> get_vertices() const; //Returns a vector containing all the vertices.
	vector<vertex<T>> get_neighbours(const int&) const; //Returns a vector containing all the vertices reachable from the given vertex. The vertex is not considered a neighbour of itself.
	vector<vertex<T>> get_second_order_neighbours(const int&); // Returns a vector containing all the second_order_neighbours (i.e., neighbours of neighbours) of the given vertex.
															  // A vector cannot be considered a second_order_neighbour of itself.
	bool reachable(const int&, const int&) const; //Returns true if the second vertex is reachable from the first (can you follow a path of out-edges to get from the first to the second?). Returns false otherwise.
	bool contain_cycles() const; // Return true if the graph contains cycles (there is a path from any vertices directly/indirectly to itself), false otherwise.
	
	vector<vertex<T>> depth_first(const int&); //Returns the vertices of the graph in the order they are visited in by a depth-first traversal starting at the given vertex.
	vector<vertex<T>> breadth_first(const int&); //Returns the vertices of the graph in the order they are visisted in by a breadth-first traversal starting at the given vertex.

	directed_graph<T> out_tree(const int&) const; //Returns a tree starting at the given vertex using the out-edges. This means every vertex in the tree is reachable from the root.

	vector<vertex<T>> pre_order_traversal(const int&, directed_graph<T>&); // returns the vertices in the visiting order of a pre-order traversal of the tree starting at the given vertex.
	vector<vertex<T>> in_order_traversal(const int&, directed_graph<T>&); // returns the vertices in the visiting order of an in-order traversal of the tree starting at the given vertex.
	vector<vertex<T>> post_order_traversal(const int&, directed_graph<T>&); // returns the vertices in ther visitig order of a post-order traversal of the tree starting at the given vertex.

	vector<vertex<T>> significance_sorting(); // Return a vector containing a sorted list of the vertices in descending order of their significance.

};

// Define all your methods down here (or move them up into the header, but be careful you don't double up). If you want to move this into another file, you can, but you should #include the file here.
// Although these are just the same names copied from above, you may find a few more clues in the full method headers.
// Note also that C++ is sensitive to the order you declare and define things in - you have to have it available before you use it.

template <typename T>
directed_graph<T>::directed_graph() 
{
	//intial capacity for matrix (100*100)
	int initial_cap = 100;
	//Intialise the adjacency matrix
	adj_matrix.resize(initial_cap);
	//Makes a square matrix using double "for loop"
	for(int i=0; i<adj_matrix.size(); i++)
	{
		adj_matrix[i].resize(initial_cap);
		for(int j=0; j<<adj_matrix[i].size(); j++)
		{
			adj_matrix[i][j] = 0; //when the edge weight is 0 it means there is no realtionship between vertex [i] and [j]
		}
	}

	vertices_weights.resize(initial_cap);
	for(int i=0; i<vertices_weights.size() ;i++)
	{
		vertices_weights[i] = 0; //when the weight is 0 it means the vertex is not in the list
	}
}

template <typename T>
directed_graph<T>::~directed_graph() {}

template <typename T>
void directed_graph<T>::increaseCapacity()
{
	//make new capacity *2 as old one
	int old_cap = adj_matrix.size();
	int new_cap = old_cap * 2;

	//resize the adjacent matrix
	adj_matrix.resize(new_cap);
	for(int i=0; i<adj_matrix.size(); i++)
	{
		adj_matrix[i].resize(new_cap);
	}
	//Initialise the adjacency matrix
	for(int i=0; i<old_cap; i++)
	{
		for(int j=old_cap; j<new_cap; j++)
		{
			adj_matrix[i][j] = 0; //Intialise all edges as 0 as there are no connections yet
		}
	}

	//Expand the old vertext weights into new ones
	vertices_weights.resize(new_cap);
	//Starts with the last capacity
	for(int i=old_cap; i<new_cap ;i++)
	{
		vertices_weights[i] = 0; //when the weight is 0 it means the vertex is not in the list
	}
}

template <typename T>
bool directed_graph<T>::contains(const int& u_id) const 
{
	//Checks the vertice_weight vector with the index of u_id
	if(vertices_weights[u_id] > 0)
	{
		//All vertex needs to have weight, if it does not have a weight (0) than it does not exist
		return true;
	}
	//returns false when vertex pointed has no weight (vertex does not exist)
	return false; 
}

template <typename T>
bool directed_graph<T>::adjacent(const int& u_id, const int& v_id) const 
{
	//Checks if the 2 vectors
	if(adj_matrix[u_id][v_id] > 0)
	{
		return true;
	} 
	return false; 
}

template <typename T>
void directed_graph<T>::add_vertex(const vertex<T>& u) 
{
	//when the id is larger than the last index of vector increase the vertex_weights and adj_matrix
	while(u.id > adj_matrix.size()-1)
	{
		increaseCapacity();
	}
	//vertices_weight[index = u.id] = (stores the weight of the vertex)
	vertices_weights[u.id] = u.weight;
}

template <typename T>
void directed_graph<T>::add_edge(const int& u_id, const int& v_id, const T& weight) 
{
	if(contains(u_id) && contains(v_id)) //Checks if the vertices id (u_id & v_id) exist in the first place
	{
		if(weight > 0) //The weight added must be 0 or because all edges have weight in a "directed weighted graph"
		{
			adj_matrix[u_id][v_id] = weight; //Adds the weight value to the adjacency matrix coordinate/position Adj_Matrix[u_id][v_id]
		}
	}
}

template <typename T>
void directed_graph<T>::remove_vertex(const int& u_id) 
{
	//Checks of the vertex exist in the first place
	if(contains(u_id))
	{
		//Changes the value of the vertex in the vertices_weight to 0
		//A vertex with 0 weight does not exist 
		vertices_weights[u_id] = 0;
	}
}

template <typename T>
void directed_graph<T>::remove_edge(const int& u_id, const int& v_id) 
{
	//Checks if the vertex u_id and v_id exist
	if(contains(u_id) && contains(v_id))
	{
		//Changes the value of the edge at coordinate/position [u_id][v_id] to 0
		//0 means there are no connections between u_id and v_id
		adj_matrix[u_id][v_id] = 0;
	}
}

template <typename T>
size_t directed_graph<T>::in_degree(const int& u_id) const 
{
	//variable to store the number of in_degree (initialised as 0)
	int num_of_in_degree = 0;
	//Iterates through the adjacency matrix at a specific row to find any edge value more than 0
	//Directly connected to the vertex is adjacent
	for(int i=0; i < adj_matrix.size(); i++)
	{
		//iterating through adj_matrix[..][u_id] to find a value more than 0 (means incoming to vertex u_id)
		//adj_matrix[i][u_id] = adj_matrix["From"]["To"] =adj_matrix["Iterates"]["Remains the same"]
		if(adj_matrix[i][u_id] > 0)
		{
			//Everytime an incoming degree is found increment the num_of_degree
			//Eventually adds up to the total of in_degree for the u_id vertex
			num_of_in_degree++;
		}
	}
	//returns the total number of in degree
	return num_of_in_degree;
}

template <typename T>
size_t directed_graph<T>::out_degree(const int& u_id) const 
{
	//variable to store the number of out_degree (initialised as 0)
	int num_of_out_degree = 0;
	//iterating through adj_matrix[u_id][..] to find a value more than 0 (means incoming to vertex u_id)
	//adj_matrix[u_id][i] = adj_matrix["From"]["To"] =adj_matrix["Remains the same"][""Iterates""]
	for(int i=0; i < adj_matrix.size(); i++)
	{
		//Everytime an out degree is found increment the num_of_out_degree
		//Eventually adds up to the total of out_degree for the u_id vertex
		if(adj_matrix[u_id][i] > 0)
		{
			//Everytime an out degree is found increment the num_of_out_degree
			//Eventually adds up to the total of in_degree for the u_id vertex
			num_of_out_degree++;
		}
	}
	//returns the total number of out degree
	return num_of_out_degree;
}

template <typename T>
size_t directed_graph<T>::degree(const int& u_id) const 
{
	//Returns the number of all degrees (out and in)
	//Calls the in_degree and out_degree function
	//Both in_degree() and out_degree() returns the total number of degrees for vertex u_id
	//Sum both of the value (in_degree and out_degree) together to get the total degree
	int num_of_degree = out_degree(u_id) + in_degree(u_id);
	//Returns the total number of degree for vertex u_id
	return num_of_degree;
}

template <typename T>
size_t directed_graph<T>::num_vertices() const 
{
	//variable to store the number of vertices (initialised as 0)
	int num_vertices = 0;
	//Iterates through the vertices_weights that keeps (the id of the vertices and the weight)
	for(int i =0; i < vertices_weights.size(); i++)
	{
		//Checks if the index "i" in vertices_weights has weight using contains() function
		if(contains(i))
		{
			//Everytime there is vertex increment the num_vertices
			//Eventually totals up to the total number of vertices
			num_vertices++;
		}
	}
	//returns the total number of vertices
	return num_vertices;
}

template <typename T>
size_t directed_graph<T>::num_edges() const 
{
	//variable to store the number of vertices (initialised as 0)
	int num_edges = 0;
	//Iterates through the adjacency matrix that keeps the edges
	for(int i=0; i<adj_matrix.size(); i++)
	{
		for(int j=0; j<adj_matrix.size(); j++)
		//Check if there is an edge in adjacency matrix coordinate i,j
		if(adj_matrix[i][j] > 0)
		{
			//Everytime there is vertex increment the num_edges
			//Eventually totals up to the total number of edges
			num_edges++;
		}
	}
	//returns the total number of edges
	return num_edges; 
}

template <typename T>
vector<vertex<T>> directed_graph<T>::get_vertices() const
{
	//Vector to store vertices that exist
	vector<vertex<T>> vertice_list;
	//Iterates through vertices_weight to find the vertices
	for(int i=0; i<vertices_weights.size(); i++)
	{
		//Checks if the index "i" in vertices_weights has weight using contains() function
		if(contains(i))
		{
			//Everytime their is a vertex
			//push the a vertext object (make vertex object at with id "i" and the weight of vertices_weight[i]
			//Basically just copying value from vertices_weights to vertices_list to return a list of vertices
			vertice_list.push_back(vertex<T>(i, vertices_weights[i]));
		}
	}
	
	//Returns vertice_list, a vector that contains all the vertices
	return vertice_list;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::get_neighbours(const int& u_id) const
{
	//Vectors to keep a list of the first order neighbour
	vector<vertex<T>> neighbour_list;
	//Checks if u_id is a vertex that exist
	if(contains(u_id))
	{
		//Iterates through the adjacency matrix to find a vertex directly connected to it
		for(int i=0; i<adj_matrix[u_id].size(); i++)
		{
			//Iterates through the adjacency matrix at index [u_id] to find a direct connection, a first order neighbour
			//Check if the edge at adjacency matrix at coordinate [u_id],[i] has an edge weight of more than 0
			//Edge weight == 0 neans that there is no edge
			if(adj_matrix[u_id][i] > 0)
			{
				//Everytime there is an edge push the vertice to the neighbour list
				neighbour_list.push_back(vertex<T>(i, vertices_weights[i]));
			}
		}
	}
	//Returns the neigbour_list, a vector that contains all vertex that are connected directly
	return neighbour_list; 
}

template <typename T>
vector<vertex<T>> directed_graph<T>::get_second_order_neighbours(const int& u_id) 
{
	//a vector, firs_order to keep the value of the first order neigbours
	//Stores the results of the get_neigbours(u_id) that gets the first order neigbours of u_id
	vector<vertex<T>> first_order = get_neighbours(u_id);
	//A vector, second_order to keep the value of the second order neigbour
	vector<vertex<T>> second_order;
	//A vector, in_second_order to keep ensure no duplicates exist in the second_order
	//This is because a vertex may be a second_order neighbour to 2 different vertex
	vector<bool>in_second_order (adj_matrix.size(), false);

	//for each vertex in first_order (first order neighbour) get their neighbour (the second order neigbour of u_id)
	//keyword auto to determine the datatype automatically
	for(auto vertex_first :first_order)
	{
		//Iterates through the adjacent matrix looking for the negihbour of the first order neighbour (second order neigbour of u_id)
		for(int i=0; i<adj_matrix[u_id].size(); i++)
		{
			//Checks if the their is an edge between the firs_order neighbour vertex and "i"
			//Makes sure that "i" is not equal to u_id as a vertex cannot be it's own neigbour
			//Checks if the vertex is not in the second_order already, this is to prevent duplicates
			if(adj_matrix[vertex_first.id][i] > 0 && i != u_id && !(in_second_order[i]))
			{
				//Everytime there is a "new" second neighbour push the vertex to the list
				//The value will be the id of the vertex "i" and the weight from the vertice_weights at index "i"
				second_order.push_back(vertex<T>(i, vertices_weights[i]));
				//Mark the vertex as true, because it has been added to the second_order, this is to prevent duplicates
				//Because a vertex can be a second neighbour to more than 1 vertex thus, duplicates might come up
				in_second_order[i] = true;
			}
		}
	}
	return second_order;
}

template <typename T>
bool directed_graph<T>::reachable(const int& u_id, const int& v_id) const 
{
	//Boolean variable to return
	bool reachable = false;
	//Checks if u_id and v_id exist in the graph
	if(contains(u_id) && contains(v_id))
	{
		//vertice list to keep 
		vector<int> vertice_list;
		vector<int> q;
		vector<bool> visited(adj_matrix.size(), false);

		//push the u_id (start) to q (the exploration list)
		q.push_back(u_id);
		//Mark u_id as visisted
		visited[u_id] = true;

		int visit;
		while(!q.empty())
		{
			visit = q[0];
			vertice_list.push_back(visit);
			q.erase(q.begin());

			for(int i=0; i < adj_matrix.size(); i++)
			{
				if(adj_matrix[visit][i] > 0 && (!visited[i]))
				{
					q.push_back(i);
					visited[i] = true;
				}
			}
		}
		//iterate through the vertice list that contains all vertex visited throguh BFS start at u_id (all possible visited vertex)
		//Checks if v_id is in the list, if it's in the list than it is reachable
		for(int i=0; i<vertice_list.size(); i++)
		{
			//Checks if there is a v_id value in the vertice_list made
			if(vertice_list[i] == v_id)
			{
				//v_id is reachable if it's in vertice_list
				//Change reachable variable to true because it's reachable
				reachable = true;
			}
		}
	}
	//Returns a boolean variable wheter the v_id is reachable from u_id
	return reachable;
}


template <typename T>
bool directed_graph<T>::contain_cycles() const
{
	//Boolean variable to keep wheter the graph is cyclic
	bool isCyclic = false;
	//Get and store a list of all vertices in the graph
	vector<vertex<T>> vertex_list = get_vertices();
	vector<vertex<T>> neighbour;
	//For every vertex/node on the graph get it's respective neigbours
	//Check if there is a "back edge" the child node can reach the parent node
	//Check the neighbour ("child") of each vertex ("parent") and see if they are reacheble through it
	//If they are reachable than their is a "back edge" thus, there is a cycle in the graph; 
	for(auto vertex:vertex_list)
	{
		//Checks if the the vertex exist
		if(contains(vertex.id))
		{
			//Gets the neighbour of the vertex
			neighbour = get_neighbours(vertex.id);
			for(auto n:neighbour)
			{
				//For every neigbour of the vertex check if it is reachable
				//Detect if their is a backedge
				if(reachable(n.id, vertex.id))
				{
					//If the the neigbour can reach the vertex than there is a cycle
					isCyclic = true;
				}
			}
		}
	}
	//Returns isCyclic, the variable that true of there is a cycle and false if there is no cycle
	return isCyclic;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::depth_first(const int& u_id) 
{
	//dfs_vertice_list, a vector to keep the vertex visited during dfs
	vector<vertex<T>> dfs_vertice_list;
	//to_be_visited, a stack that contains a list of vertex that needs to be explored
	stack<int> to_be_visited;
	//visited, a boolean vector to keep track of visited vector (intialised as all false because nothing is visited yet, initially)
	//The size of the vector is the same as the adj_matrix so, it can fit all the vertex
	vector<bool>visited(adj_matrix.size(), false);
	
	//push u_id to the stack, to_be_visited to start the dfs
	to_be_visited.push(u_id);

	//Keep exploring until there is no more vertex to visit
	//Loop that keeps going until the the stack, to_be_visited is empty
	while(!to_be_visited.empty())
	{
		//varibale to keep the last value inputed into the stack, to_be_visited (top value)
		//This is because the stack uses the Last In First Out Principle
		int visit = to_be_visited.top();
		//Pop/delete the top of the stack so, we can explore the next vertex
		to_be_visited.pop();
		//Checks if the vertex has been visited before (this is to prevent an infinite loop)
		if(!visited[visit])
		{
			//Set the vertext as visited so, it won't be visited again
			visited[visit] = true;
			//Push the vertex visited to the vector, dfs_vertice_list, to keep track of vertices visited by dfs
			dfs_vertice_list.push_back(vertex<T>(visit, vertices_weights[visit]));
			//Iterate through the adjacency matrix from the root (starting from the end), exploring all adjacent vertex
			//starts from the last index of the matrix = adj_matrix.size() and move down i--
			for(int i=adj_matrix.size(); i != 0; i--)
			{
				//Check for the existence of the edge
				if(adj_matrix[visit][i-1] > 0)
				{
					//push the vertex to the top of the stack, to_be_visited so, it will be visited
					to_be_visited.push(i-1);
				}
			}
		}
	}
	//Returns a list of vertex that has been visited during dfs
	return dfs_vertice_list; 
}

template <typename T>
vector<vertex<T>> directed_graph<T>::breadth_first(const int& u_id) 
{
	//vertice list to keep the list of vertices visisted in breadth first search
	vector<vertex<T>> vertice_list;
	//Vector to keep the list of vertice that needs to be visit
	vector<int> to_be_visited;
	//Vector to mark a node as visited/not visited (true/false), this prevents loops
	vector<bool> visited(adj_matrix.size(), false);

	//push the u_id (start) to queue (the exploration list)
	to_be_visited.push_back(u_id);
	//Mark u_id (the first vertex visited) as visisted
	visited[u_id] = true;

	//Variable to store the first element of the queue
	int visit;
	//Run until their is no vertex left unexplored
	while(!to_be_visited.empty())
	{
		//visit the first vertex in the queue (because queue is first in first out)
		visit = to_be_visited[0];
		//Push vertex to the vertice_list after the vertex is visited
		vertice_list.push_back(vertex<T>(visit, vertices_weights[visit]));
		//pop the first element in the queue so, in the next iteration we visit a another vertex
		to_be_visited.erase(to_be_visited.begin());

		//Iterates through the adjacency matrix to find all adjacent vertex to the vertex taken by the queue
		for(int i=0; i < adj_matrix.size(); i++)
		{
			//Add the vertex only to the queue if the vertex has never been visited and has a connection to the current vertex 
			if(adj_matrix[visit][i] > 0 && (!visited[i]))
			{
				//Add the adjacent vertex to the queue of vertex to_be_visited
				//vertex are keep being added to the queue until there is no more vertex that can be explored or reached
				to_be_visited.push_back(i);
				//Mark the adjacent vertex as visited to prevent an infinite loop
				visited[i] = true;
			}

		}
	}
	return vertice_list;
}

template <typename T>
directed_graph<T> directed_graph<T>::out_tree(const int& u_id)  const
{
	//Create a directed_graph called tree to store the out_tree
	directed_graph<T> tree;
	//Queue type integer to store the vertex that needs to be explored next
	queue<int> to_be_visited;
	//Insert u_id into the queue to start the queue
	to_be_visited.push(u_id);
	//Insert u_id to the tree as the first node (root)
	//Add by first creating the vertex and adding the u_id and the weight of the node (from the main graph)
	vertex<int> root(u_id, vertices_weights[u_id]);
	//Add the root to the tree by using the "add_vertex" method (part of the directed_graph class)
	tree.add_vertex(root);

	//Keep exploring until there is no more node to explore on the queue
	while(!to_be_visited.empty())
	{
		//variable to keep the first element of the queue
		//element = node that is going to be visit
		int visit = to_be_visited.front();
		//Remove the element from the queue (so, the next node can be explored in the next iteration)
		to_be_visited.pop();

		//Explore all the adjacent node to the current node
		//Check if all adjacent nodes are in the tree
		for(int i=0; i<adj_matrix[i].size(); i++)
		{
			if(adj_matrix[visit][i] > 0)
			{
				//Checks if the node exist already in the tree
				//Using contains method to check if the node is already included in the tree
				if(!tree.contains(i))
				{
					//Create the node/vertex that is going to be added
					//Add the vertice_id and the weight based on vertice_weights (from main graph)
					vertex<int> v(i, vertices_weights[i]);
					//Add the node/vertex to the tree using the add_vertex method (method from directed_graph class)
					tree.add_vertex(v);
					//Connect the adjacent node with the currently visited node
					//Connect by adding edges using the add_edge method (method from directed_graph class)
					//Pass the "visit" and "i" as the nodes and get the weight of the node from the adj_matrix of the main graph
					tree.add_edge(visit, i, adj_matrix[visit][i]);
					//Push the adjacent node to the queue for furthe exploration
					to_be_visited.push(i);
				}
			}
		}
	}
	//returns a directed_graph that contains the out_tree created
	return tree;
}

template <typename T>
vector<vertex<T>> directed_graph<T>::pre_order_traversal(const int& u_id, directed_graph<T>& mst) 
{
	//dfs_vertice_list, a vector to keep the vertex visited during dfs
	vector<vertex<T>> dfs_vertice_list;
	//to_be_visited, a stack that contains a list of vertex that needs to be explored
	stack<int> to_be_visited;
	//visited, a boolean vector to keep track of visited vector (intialised as all false because nothing is visited yet, initially)
	//The size of the vector is the same as the adj_matrix so, it can fit all the vertex
	vector<bool>visited(mst.adj_matrix.size(), false);
	
	//push u_id to the stack, to_be_visited to start the dfs
	to_be_visited.push(u_id);

	//Keep exploring until there is no more vertex to visit
	//Loop that keeps going until the the stack, to_be_visited is empty
	while(!to_be_visited.empty())
	{
		//varibale to keep the last value inputed into the stack, to_be_visited (top value)
		//This is because the stack uses the Last In First Out Principle
		int visit = to_be_visited.top();
		//Pop/delete the top of the stack so, we can explore the next vertex
		to_be_visited.pop();
		//Checks if the vertex has been visited before (this is to prevent an infinite loop)
		if(!visited[visit])
		{
			//Set the vertext as visited so, it won't be visited again
			visited[visit] = true;
			//Push the vertex visited to the vector, dfs_vertice_list, to keep track of vertices visited by dfs
			dfs_vertice_list.push_back(vertex<T>(visit, vertices_weights[visit]));
			//Iterate through the adjacency matrix from the root (starting from the end), exploring all adjacent vertex
			//starts from the last index of the matrix = adj_matrix.size() and move down i--
			for(int i=mst.adj_matrix.size(); i != 0; i--)
			{
				//Check for the existence of the edge
				if(mst.adj_matrix[visit][i-1] > 0)
				{
					//push the vertex to the top of the stack, to_be_visited so, it will be visited
					to_be_visited.push(i-1);
				}
			}
		}
	}
	//Returns a list of vertex that has been visited during dfs
	return dfs_vertice_list; 
}

template <typename T>
vector<vertex<T>> directed_graph<T>::in_order_traversal(const int& u_id, directed_graph<T>& mst) 
{
	vector<vertex<T>> vertices_list;
	return vector<vertex<T>>(); 
}

template <typename T>
vector<vertex<T>> directed_graph<T>::post_order_traversal(const int& u_id, directed_graph<T>& mst) { return vector<vertex<T>>(); }

template <typename T>
vector<vertex<T>> directed_graph<T>::significance_sorting() 
{
	//Gets all vertices along with their weight and id using get.vertices() and stores them in the vertices_sorted
	vector<vertex<T>> vertices_sorted = get_vertices();
	//use sort from the C++ STL (Algorithm) to sort the vertices based on their weight
	sort(vertices_sorted.begin(), vertices_sorted.end(), [](auto& first_element, auto& second_element) {return first_element.weight > second_element.weight;});
	//returns the vertices sorted
	return vertices_sorted;
}

#endif