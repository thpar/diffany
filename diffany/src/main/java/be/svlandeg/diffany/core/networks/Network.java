package be.svlandeg.diffany.core.networks;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.algorithms.NetworkCleaning;
import be.svlandeg.diffany.core.io.NetworkIO;

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
	protected Set<String> nodeAttributes;

	protected int ID;
	protected String name;

	/**
	 * Create a new network with a specific name, ID, and sets of nodes and edges.
	 * The ID should be unique within a project, the name is preferably unique as well.
	 * All source and target nodes of each edge will be automatically added to the internal set of nodes.
	 * 
	 * @param name the name of this network 
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 * @param nodes the nodes of this network (should all contain the correct attributes if there are any defined!)
	 * @param edges the edges of this network
	 */
	public Network(String name, int ID, Set<String> nodeAttributes, Set<Node> nodes, Set<Edge> edges)
	{
		this.name = name;
		this.ID = ID;
		if (nodeAttributes != null)
		{
			this.nodeAttributes = nodeAttributes;
		}
		else
		{
			this.nodeAttributes = new HashSet<String>();
		}
		setNodesAndEdges(nodes, edges);
	}

	/**
	 * Create a new network with an empty set of nodes and edges.
	 * 
	 * @param name the name of this network 
	 * @param ID the unique identifier of this network (should be enforced to be unique within one project)
	 * @param nodeAttributes the required node attribute names for this network - can be left empty or null
	 */
	public Network(String name, int ID, Set<String> nodeAttributes)
	{
		this(name, ID, nodeAttributes, new HashSet<Node>(), new HashSet<Edge>());
	}

	/**
	 * Return the name of this network
	 * @return the name of this network
	 */
	public String getName()
	{
		return name;
	}
	
	/**
	 * Return the ID of this network (which should be unique within one project)
	 * @return the ID of this network
	 */
	public int getID()
	{
		return ID;
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
	 * Get the set of edges in this network.
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
				if (!e.isSymmetrical())
				{
					if (e.getSource().ID.equals(source.ID) && e.getTarget().ID.equals(target.ID))
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
		String sourceID = source.ID;
		String targetID = target.ID;
		
		return getAllEdges(sourceID, targetID);
	}
	
	/**
	 * Get all edges (both symmetric and assymetric) in this network between two specific nodes. 
	 * In case there are symmetrical edges in this network between target-source, these will be added too.
	 * 
	 * @param sourceID the required source node, defined by ID
	 * @param targetID the required target node, defined by ID
	 * @return the set of edges between these two nodes (can be empty, but not null)
	 */
	public Set<Edge> getAllEdges(String sourceID, String targetID)
	{
		Set<Edge> resultEdges = new HashSet<Edge>();
		if (sourceID != null && targetID != null)
		{
			for (Edge e : edges)
			{
				String edgeSourceID = e.getSource().ID;
				String edgeTargetID = e.getTarget().ID;
				
				if (edgeSourceID.equals(sourceID) && edgeTargetID.equals(targetID))
				{
					resultEdges.add(e);
				}
				else if (e.isSymmetrical() && edgeSourceID.equals(targetID) && edgeTargetID.equals(sourceID))
				{
					resultEdges.add(e);
				}
			}
		}
		return resultEdges;
	}

	/**
	 * Get all edge definitions in this network between two specific nodes. 
	 * In case symmetry=true and there are symmetrical edges in this network between target-source, these will be added too.
	 * 
	 * @param source the required source node 
	 * @param target the required target node
	 * @param symmetry the symmetrical state of the edges, defining the result set to either be directed or symmetrical
	 * @return the set of edge definitions between these two nodes (can be empty, but not null)
	 */
	public Set<EdgeDefinition> getAllEdgeDefinitions(Node source, Node target, boolean symmetry)
	{
		Set<EdgeDefinition> resultEdgeDefinitions = new HashSet<EdgeDefinition>();
		Set<Edge> resultEdges = getAllEdges(source, target);
		for (Edge e : resultEdges)
		{
			if (e.isSymmetrical() == symmetry)
			{
				resultEdgeDefinitions.add(e.def);
			}
		}
		return resultEdgeDefinitions;
	}

	
	
	/**
	 * This method removes all unconnected nodes from the network - can not be reverted!
	 */
	public void removeUnconnectedNodes()
	{
		setNodesAndEdges(edges);
	}
	
	/**
	 * Define a (new) set of edges for this network, overwriting previous data.
	 * The node set will be (only) those appearing as source or target in the given edges.
	 * 
	 * @param edges the edges of this network (implicitly also defining the nodes)
	 */
	public void setNodesAndEdges(Set<Edge> edges)
	{
		setNodesAndEdges(new HashSet<Node>(), edges);
	}

	/**
	 * Define a (new) set of nodes for this network, overwriting previous data.
	 * All source and target nodes of each edge will be automatically added to the internal set of nodes.
	 * 
	 * @param edges the edges of this network
	 * @param nodes the nodes of this network
	 */
	public void setNodesAndEdges(Set<Node> nodes, Set<Edge> edges)
	{
		this.nodes = new HashSet<Node>();
		for (Node n : nodes)
		{
			addNode(n);
		}
		this.edges = new HashSet<Edge>();
		for (Edge e : edges)
		{
			addEdge(e);
		}
	}

	/**
	 * Add an edge to this network, automatically also adding its source and target nodes if needed.
	 * If the edge was already present, nothing happens. An additional cleaning step (with {@link NetworkCleaning}) should be performed to remove redundant edges. 
	 * 
	 * @param edge a new adge in this network
	 */
	public void addEdge(Edge edge)
	{
		edges.add(edge);

		addNode(edge.getSource());
		addNode(edge.getTarget());
	}

	/**
	 * Remove an edge from this network, leaving its source and target nodes otherwise untouched (and perhaps unconnected).
	 * @param edge the edge that should be removed from the network
	 */
	public void removeEdge(Edge edge)
	{
		edges.remove(edge);
	}
	
	/**
	 * Get all node attributes required in this network.
	 * 
	 * @return the set of all node attributes
	 * @throws IllegalArgumentException when the current network already contains some nodes
	 */
	public Set<String> getAllNodeAttributes()
	{
		return nodeAttributes;
	}
	
	/**
	 * Add a new required node attribute to this network. This is only possible when the current network is empty (i.e. there are no nodes or edges).
	 * 
	 * @param nodeAttribute a new node attribute
	 * @throws IllegalArgumentException when the current network already contains some nodes
	 */
	public void addNodeAttribute(String nodeAttribute)
	{
		if (nodes.size() > 0)
		{
			String errormsg = "Can not add node attributes when the network already contains nodes!";
			throw new IllegalArgumentException(errormsg);
		}
		nodeAttributes.add(nodeAttribute);
	}
	
	/**
	 * Define the node attributes by taking only those that are in common from a given set of networks.
	 * This is used for instance at the creation of differential/consensus networks.
	 * 
	 * @param originalNetworks the original networks which are scanned for their node attributes
	 */
	protected void defineCommonAttributes(Set<Network> originalNetworks)
	{
		Map<String, Integer> countedAttributes = new HashMap<String, Integer>();
		for (Network n : originalNetworks)
		{
			for (String att : n.getAllNodeAttributes())
			{
				if (!countedAttributes.containsKey(att))
				{
					countedAttributes.put(att, 0);
				}
				countedAttributes.put(att, countedAttributes.get(att) + 1);
			}
		}
		for (String att : countedAttributes.keySet())
		{
			int count = countedAttributes.get(att);
			if (count == originalNetworks.size())
			{
				addNodeAttribute(att);
			}
		}
		
	}

	/**
	 * Add a new (unconnected) node to this network. If it was already present, nothing happens. The comparison is done with the internal NodeMapper object.
	 * The comparison is made through the NodeMapper object of this Network.
	 * 
	 * @param node a new node in this network
	 * @throws IllegalArgumentException when the given node does not contain the required attributes, as specified by this network
	 */
	public void addNode(Node node)
	{
		boolean contained = false;
		for (Node n : nodes)
		{
			if (n.ID.equals(node.ID))
			{
				contained = true;
			}
		}
		if (!contained)
		{
			nodes.add(node);
		}
		for (String nodeAttribute : nodeAttributes)
		{
			Object value = node.getAttribute(nodeAttribute);
			if (value == null)
			{
				String errormsg = "The node " + node.ID + " does not contain the required attribute " + nodeAttribute + "!";
				throw new IllegalArgumentException(errormsg);
			}
		}
	}
	
}
