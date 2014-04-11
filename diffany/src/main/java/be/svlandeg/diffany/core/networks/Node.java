package be.svlandeg.diffany.core.networks;

/**
 * Class that represents a node in a network. A node can be virtual, i.e. non-existing, for the purpose of modeling
 * unknown interaction partners which may still convey important information.
 * 
 * @author Sofie Van Landeghem
 */
public class Node
{
	
	protected String name;
	protected boolean virtual;
	
	/**
	 * Create a new (non-virtual) node with a specific name
	 * @param name the name of this node - should be unique within a network!
	 */
	public Node(String name) throws IllegalArgumentException
	{
		if (name == null)
		{
			String errormsg = "The name of a node should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.name = name;
		virtual = false;
	}
	
	/**
	 * Create a new node with a specific name and which is virtual or not.
	 * @param name the name of this node - should be unique within a network!
	 * @param virtual whether this node is virtual
	 */
	public Node(String name, boolean virtual) throws IllegalArgumentException
	{
		if (name == null)
		{
			String errormsg = "The name of a node should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.name = name;
		this.virtual = virtual;
	}
	
	/**
	 * Return the name of this node, in its original (unnormalized) form.
	 * @return the name of this node
	 */
	public String getName()
	{
		return getName(false);
	}

	/**
	 * Return the name of this node, either normalized or in original form.
	 * 
	 * @param normalized whether or not the name should be converted to lower case.
	 * @return the name of this node
	 */
	public String getName(boolean normalized)
	{
		if (normalized)
		{
			return name.toLowerCase();
		}
		return name;
	}
	
	/**
	 * Check whether or not this node is virtual or not.
	 * E.g. visualisations will want to hide virtual nodes.
	 * 
	 * @return whether or not this node is virtual
	 */
	public boolean isVirtual()
	{
		return virtual;
	}
	
	@Override
	public String toString()
	{
		return "Node " + name;
	}
}
