import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class Graph {

  // Keep a fast index to nodes in the map
  private Map<Integer, Vertex> vertices;

  /**
   * Construct an empty Graph with a map. The map's key is the name of a vertex
   * and the map's value is the vertex object.
   */
  public Graph() {
    vertices = new HashMap<>();
  }

  /**
   * Adds a vertex to the graph. Throws IllegalArgumentException if two vertices
   * with the same name are added.
   * 
   * @param v
   *          (Vertex) vertex to be added to the graph
   */
  public void addVertex(Vertex v) {
    if (vertices.containsKey(v.name))
      throw new IllegalArgumentException("Cannot create new vertex with existing name.");
    vertices.put(v.name, v);
  }

  /**
   * Gets a collection of all the vertices in the graph
   * 
   * @return (Collection<Vertex>) collection of all the vertices in the graph
   */
  public Collection<Vertex> getVertices() {
    return vertices.values();
  }

  /**
   * Gets the vertex object with the given name
   * 
   * @param name
   *          (String) name of the vertex object requested
   * @return (Vertex) vertex object associated with the name
   */
  public Vertex getVertex(int name) {
    return vertices.get(name);
  }

  /**
   * Adds a directed edge from vertex u to vertex v
   * 
   * @param nameU
   *          (String) name of vertex u
   * @param nameV
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addEdge(int nameU, int nameV, Double cost) {
    if (!vertices.containsKey(nameU))
      throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
    if (!vertices.containsKey(nameV))
      throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
    Vertex sourceVertex = vertices.get(nameU);
    Vertex targetVertex = vertices.get(nameV);
    Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
    sourceVertex.addEdge(newEdge);
  }

  /**
   * Adds an undirected edge between vertex u and vertex v by adding a directed
   * edge from u to v, then a directed edge from v to u
   * 
   * @param name
   *          (String) name of vertex u
   * @param name2
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addUndirectedEdge(int name, int name2, double cost) {
    addEdge(name, name2, cost);
    addEdge(name2, name, cost);
  }


  /**
   * Computes the euclidean distance between two points as described by their
   * coordinates
   * 
   * @param ux
   *          (double) x coordinate of point u
   * @param uy
   *          (double) y coordinate of point u
   * @param vx
   *          (double) x coordinate of point v
   * @param vy
   *          (double) y coordinate of point v
   * @return (double) distance between the two points
   */
  public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
    return Math.sqrt(Math.pow(ux - vx, 2) + Math.pow(uy - vy, 2));
  }

  /**
   * Computes euclidean distance between two vertices as described by their
   * coordinates
   * 
   * @param u
   *          (Vertex) vertex u
   * @param v
   *          (Vertex) vertex v
   * @return (double) distance between two vertices
   */
  public double computeEuclideanDistance(Vertex u, Vertex v) {
    return computeEuclideanDistance(u.x, u.y, v.x, v.y);
  }

  /**
   * Calculates the euclidean distance for all edges in the map using the
   * computeEuclideanCost method.
   */
  public void computeAllEuclideanDistances() {
    for (Vertex u : getVertices())
      for (Edge uv : u.adjacentEdges) {
        Vertex v = uv.target;
        uv.distance = computeEuclideanDistance(u.x, u.y, v.x, v.y);
      }
  }

  Random rd = new Random();
  public int numberOfVertices = 0;

  public void generateRandomVertices(int n) {
    numberOfVertices = n;
    vertices = new HashMap<>(); // reset the vertex hashmap
    for(int i = 0; i < n; i++) {
      Vertex newVertex = new Vertex(i, rd.nextInt(100), rd.nextInt(100));
      addVertex(newVertex);
      for(int j = i; j>=0 ; j--) {
        addUndirectedEdge(i, j, 1.00);
      }
    }
    computeAllEuclideanDistances(); // compute distances
  }

  public List<Edge> nearestNeighborTsp() {

    Vertex root = this.getVertex(rd.nextInt(numberOfVertices));
    List<Edge> nearestEdges = new ArrayList<>();
    Vertex tmp = root;
    tmp.known = true;

    while(!isACompleteCycle()) {
      double minimumCost = Double.POSITIVE_INFINITY;
      Edge lowestCostEdge = null;
      for (Edge edge: tmp.adjacentEdges) {
        if(edge.target.known == false) {
          if(edge.distance < minimumCost) {
            lowestCostEdge = edge;
            minimumCost = edge.distance;
          }
        }
      }
      nearestEdges.add(lowestCostEdge);
      tmp = lowestCostEdge.target;
      tmp.known = true;
    }
    for(Edge e: root.adjacentEdges) {
      if(e.target == tmp) {
         nearestEdges.add(e);
      }
    }
    return nearestEdges;
  }

  public boolean isACompleteCycle() {
    for(Vertex v: this.getVertices()) {
      if(v.known == false) {
        return false;
      }
    }
    return true;
  }

  public List<List<Vertex>> getPermutations(List<Vertex> list) {
    List<List<Vertex>> result = new ArrayList<>();
    result.add(new ArrayList<>());

    for (int i = 0; i < list.size(); i++) {
      List<List<Vertex>> allPermutations = new ArrayList<>();
      for (List<Vertex> permutation : result) {
        for (int j = 0; j < permutation.size() + 1; j++) {
          permutation.add(j, list.get(i));
          allPermutations.add(new ArrayList<>(permutation));
          permutation.remove(j);
        }
      }
      result = new ArrayList<>(allPermutations);
    }
    return result;
  }

  public List<Edge> bruteForceTsp() {
    double lowestCost = Double.POSITIVE_INFINITY;
    List<Edge> best = new ArrayList<>();
    for(List<Vertex> singlePath : getPermutations(new ArrayList<>(this.getVertices()))) {
      double pathCost = 0;
      List<Edge> path = new ArrayList<>();
      for(int i = singlePath.size()-1; i>=0; i--) {
        for(Edge edge: singlePath.get(i).adjacentEdges) {
          if(i==0 && edge.target == singlePath.get(singlePath.size()-1) ) {
            pathCost += edge.distance;
            path.add(edge);
          }
          else if(i!=0 && edge.target == singlePath.get(i-1)) {
            pathCost += edge.distance;
            path.add(edge);
          }
        }
      }
      if(pathCost<lowestCost) {
        best = path;
        lowestCost = pathCost;
      }
    }
    return best;
  }

  /**
   * Prints out the adjacency list of the graph for debugging
   */
  public void printAdjacencyList() {
    for (int u : vertices.keySet()) {
      StringBuilder sb = new StringBuilder();
      sb.append(u);
      sb.append(" -> [ ");
      for (Edge e : vertices.get(u).adjacentEdges) {
        sb.append(e.target.name);
        sb.append("(");
        sb.append(e.distance);
        sb.append(") ");
      }
      sb.append("]");
      System.out.println(sb.toString());
    }
  }
}
