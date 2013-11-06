package be.svlandeg.cytoscape.diffany.concepts;

/**
 * Abstract class that represents a node in a network.
 * @author Sofie Van Landeghem
 *
 */
public abstract class Node
{
	
	protected String name;
	
	/**
	 * Create a new node with a specific name
	 * @param name the name of this node
	 */
	public Node(String name)
	{
		this.name = name;
	}
	
	/**
	 * Return the name of this node
	 * @return the name of this node
	 */
	public String getName()
	{
		return name;
	}

}
