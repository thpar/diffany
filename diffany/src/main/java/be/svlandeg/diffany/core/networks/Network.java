package be.svlandeg.diffany.core.networks;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.io.NetworkIO;
import be.svlandeg.diffany.core.semantics.NodeMapper;

/**
 * Abstract class that represents a network: a collection of edges and nodes
 * All source and target nodes of the edges are present in the collection of nodes,
 * but not all nodes have to be connected with edges.
 * 
 * Network data can be saved and loaded through the {@link NetworkIO} class.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class Network
{

	protected Set<Node> nodes; // ensure this set is kept consistent with the edge set!
	protected Set<Edge> edges;

	protected String name;
	protected NodeMapper nm;
	

	/**
	 * Create a new network with a specific name and sets of nodes and edges.
	 * All source and target nodes of each edge will be automatically added to the internal set of nodes.
	 * @param name the name of this network (should be enforced to be unique within one project)
	 * @param nodes the nodes of this network
	 * @param edges the edges of this network
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 */
	public Network(String name, Set<Node> nodes, Set<Edge> edges, NodeMapper nm)
	{
		if (nm == null)
		{
			String errormsg = "Please define a proper NodeMapper object!";
			throw new IllegalArgumentException(errormsg);
		}
		
		this.name = name;
		this.nm = nm;
		setNodesAndEdges(nodes, edges);
	}

	
	/**
	 * Create a new network with an empty set of nodes and edges.
	 * @param name the name of this network (should be enforced to be unique within one project)
	 * @param nm the {@link NodeMapper} object that defines equality between nodes for comparison purposes
	 */
	public Network(String name, NodeMapper nm)
	{
		this(name, new HashSet<Node>(), new HashSet<Edge>(), nm);
	}

	/**
	 * Return the name of this network (which should be unique within one project)
	 * @return the name of this network
	 */
	public String getName()
	{
		return name;
	}

	/**
	 * Obtain an easy readible string representation of this network.
	 * Ideally, the string representation includes the name of the network.
	 * @return a string representation of this network 
	 */
	public abstract String getStringRepresentation();

	/**
	 * Get the set of nodes in this network
	 * @return the set of nodes (can be empty, but not null)
	 */
	public Set<Node> getNodes()
	{
		return nodes;
	}

	/**
	 * Get the set of edges in this network
	 * @return the set of edges (can be empty, but not null)
	 */
	public Set<Edge> getEdges()
	{
		return edges;
	}

	
	/**
	 * Get all directed edges in this network between two specific nodes. 
	 * In case there are symmetrical edges in this network between source-target or target-source, these will be excluded!
	 * 
	 * @param source the required source node 
	 * @param target the required target node
	 * @return the set of edges between these two nodes (can be empty, but not null)
	 */
	public Set<Edge> getDirectedEdges(Node source, Node target)
	{
		Set<Edge> resultEdges = new HashSet<Edge>();
		if (source != null && target != null)
		{
			for (Edge e : edges)
			{
				if (! e.symmetrical)
				{
					if (nm.areEqual(e.getSource(), source) && nm.areEqual(e.getTarget(), target))
					{
						resultEdges.add(e);
					}
				}
			}
		}
		return resultEdges;
	}
	
	/**
	 * Get all edges (both symmetric and assymetric) in this network between two specific nodes. 
	 * In case there are symmetrical edges in this network between target-source, these will be added too.
	 * 
	 * @param source the required source node 
	 * @param target the required target node
	 * @return the set of edges between these two nodes (can be empty, but not null)
	 */
	public Set<Edge> getAllEdges(Node source, Node target)
	{
		Set<Edge> resultEdges = new HashSet<Edge>();
		if (source != null && target != null)
		{
			for (Edge e : edges)
			{
				if (nm.areEqual(e.getSource(), source) && nm.areEqual(e.getTarget(), target))
				{
					resultEdges.add(e);
				}
				else if (e.isSymmetrical() && nm.areEqual(e.getSource(), target) && nm.areEqual(e.getTarget(), source))
				{
					resultEdges.add(e);
				}
			}
		}
		return resultEdges;
	}
	
	/**
	 * Get all edge definitions (both symmetric and assymetric) in this network between two specific nodes. 
	 * In case there are symmetrical edges in this network between target-source, these will be added too.
	 * 
	 * @param source the required source node 
	 * @param target the required target node
	 * @return the set of edge definitions between these two nodes (can be empty, but not null)
	 */
	public Set<EdgeDefinition> getAllEdgeDefinitions(Node source, Node target)
	{
		Set<EdgeDefinition> resultEdgeDefinitions = new HashSet<EdgeDefinition>();
		Set<Edge> resultEdges = getAllEdges(source, target);
		for (Edge e : resultEdges)
		{
			resultEdgeDefinitions.add(e);
		}
		return resultEdgeDefinitions;
	}
	
	/**
	 * Get all edges in this network between two nodes with a specific name, assuming these names are unique. 
	 * In case there are symmetrical edges in this network between target-source, these will also be added if addSyms is true.
	 * 
	 * @param source the name of the required source node (normalized version in case that option is chosen)
	 * @param target the name of the required target node (normalized version in case that option is chosen)
	 * @param addSyms defined whether or not to also include symmetrical target-source edges
	 * @param normalized wehther or not the names of the nodes should be normalized before comparison to the given strings
	 * @return the set of edges between these two nodes (can be empty, but not null)
	 */
	public Set<Edge> getAllEdgesByName(String source, String target, boolean addSyms, boolean normalized)
	{
		Set<Edge> resultEdges = new HashSet<Edge>();
		if (source == null || target == null)
			return resultEdges;
		
		for (Edge e : edges)
		{
			if (e.getSource().getName(normalized).equals(source) && e.getTarget().getName(normalized).equals(target))
			{
				resultEdges.add(e);
			}
			else if (addSyms && e.isSymmetrical() && e.getSource().getName(normalized).equals(target) && e.getTarget().getName(normalized).equals(source))
			{
				resultEdges.add(e);
			}
		}
		return resultEdges;
	}

	/**
	 * Define a (new) set of nodes for this network, overwriting previous data.
	 * All source and target nodes of each edge will be automatically added to the internal set of nodes.
	 * @param edges the edges of this network
	 * @param nodes the nodes of this network
	 */
	public void setNodesAndEdges(Set<Node> nodes, Set<Edge> edges)
	{
		this.nodes = nodes;
		this.edges = new HashSet<Edge>();
		for (Edge e : edges)
		{
			addEdge(e);
		}
	}

	/**
	 * Add an edge to this network, automatically also adding its source and target nodes if needed.
	 * If the edge was already present, nothing happens 
	 * 
	 * @param edge a new adge in this network
	 */
	public void addEdge(Edge edge)
	{
		// TODO v2.0: edge comparison? Or leave it to the cleaning step?
		edges.add(edge);
		
		addNode(edge.getSource());
		addNode(edge.getTarget());
	}
	
	/**
	 * Remove an edge from this network, leaving its source and target nodes otherwise untouched.
	 * @param edge the edge that should be removed from the network
	 */
	public void removeEdge(Edge edge)
	{
		edges.remove(edge);
	}

	/**
	 * Add a new (unconnected) node to this network. If it was already present, nothing happens.
	 * The comparison is made through the NodeMapper object of this Network.
	 * 
	 * @param node a new node in this network
	 */
	public void addNode(Node node)
	{
		if (! nm.isContained(node, nodes))
		{
			nodes.add(node);
		}
	}

	


}
