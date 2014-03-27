package be.svlandeg.diffany.core.networks;

/**
 * Class that represents a node in a network.
 * 
 * @author Sofie Van Landeghem
 */
public class Node
{
	
	protected String name;
	
	/**
	 * Create a new node with a specific name
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
	
	@Override
	public String toString()
	{
		return "Node " + name;
	}

}
