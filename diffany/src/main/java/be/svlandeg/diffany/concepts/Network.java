package be.svlandeg.diffany.concepts;

import java.util.HashSet;
import java.util.Set;

/**
 * Abstract class that represents a network: a collection of edges and nodes
 * All source and target nodes of the edges are present in the collection of nodes,
 * but not all nodes have to be connected with edges.
 * @author Sofie Van Landeghem
 *
 */
public abstract class Network
{

	protected Set<Node> nodes; // ensure this set is kept consistent with the edge set!
	protected Set<Edge> edges;

	protected String name;

	/**
	 * Create a new network with a specific name and sets of nodes and edges.
	 * All source and target nodes of each edge will be automatically added to the internal set of nodes.
	 * @param name the name of this network (should be enforced to be unique within one project)
	 * @param nodes the nodes of this network
	 * @param edges the edges of this network
	 */
	public Network(String name, Set<Node> nodes, Set<Edge> edges)
	{
		this.name = name;
		setNodesAndEdges(nodes, edges);
	}

	/**
	 * Create a new network with an empty set of nodes and edges.
	 * @param name the name of this network (should be enforced to be unique within one project)
	 */
	public Network(String name)
	{
		this(name, new HashSet<Node>(), new HashSet<Edge>());
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
	 * Get all edges in this network between two specific nodes. 
	 * In case there are symmetrical edges in this network between target-source, these will also be added if addSyms is true.
	 * 
	 * @param source the required source node 
	 * @param target the required target node
	 * @param addSyms defined whether or not to also include symmetrical target-source edges
	 * @return the set of edges between these two nodes
	 */
	public Set<Edge> getAllEdges(Node source, Node target, boolean addSyms)
	{
		Set<Edge> resultEdges = new HashSet<Edge>();
		for (Edge e : edges)
		{
			if (e.getSource().equals(source) && e.getTarget().equals(target))
			{
				resultEdges.add(e);
			}
			else if (addSyms && e.isSymmetrical() && e.getSource().equals(target) && e.getTarget().equals(source))
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
	 * @param edge a new adge in this network
	 */
	public void addEdge(Edge edge)
	{
		edges.add(edge);
		nodes.add(edge.getSource());
		nodes.add(edge.getTarget());
	}

	/**
	 * Add a new (unconnected) node to this network. If it was already present, nothing happens.
	 * @param node a new node in this network
	 */
	public void addNode(Node node)
	{
		nodes.add(node);
	}

	/**
	 * Get a string representation of all edges, divided by newlines, with edges in a tabbed format.
	 * More specifically, each edge is printed as: source.name - target.name - edge.type - symmetrical - weight - negated.
	 * @return a string representation of all edges in this network, ready for printing
	 */
	public String writeEdgesTab()
	{
		String result = "";
		for (Edge e : edges)
		{
			result += e.writeToTab();
			result += System.getProperty("line.separator");
		}
		return result;
	}

	/**
	 * Remove edges in the network that are symmetrical and are represented twice (source-target and target-source).
	 * One of the two is removed only then when the type, weight and negation are all equal.
	 */
	public void removeRedundantEdges()
	{
		for (Node n1 : nodes)
		{
			String name1 = n1.getName();
			for (Node n2 : nodes)
			{
				String name2 = n2.getName();
				if (name1.compareTo(name2) < 0)
				{
					Set<Edge> edges_to = getAllEdges(n1, n2, false);
					Set<Edge> edges_back = getAllEdges(n2, n1, false);
					for (Edge et : edges_to)
					{
						for (Edge eb : edges_back)
						{
							if (!et.equals(eb))
							{
								if ((et.symmetrical && eb.symmetrical) && (et.getType().equals(eb.getType())))
								{
									if ((et.getWeight() == eb.getWeight()) && (et.isNegated() == eb.isNegated()))
									{
										edges.remove(eb);
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
