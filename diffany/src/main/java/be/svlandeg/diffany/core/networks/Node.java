package be.svlandeg.diffany.core.networks;

/**
 * Class that represents a node in a network. A node can be virtual, i.e. non-existing, for the purpose of modeling
 * unknown interaction partners which may still convey important information.
 * 
 * @author Sofie Van Landeghem
 */
public class Node
{

	protected String ID;
	protected String name;
	protected boolean virtual;
	
	/**
	 * Create a new (non-virtual) node with a specific name. The lower-case version of this name will be used as unique identifier, so ensure its unambiguity across the project!
	 * 
	 * @param name the name of this node - should be unique within a network!
	 */
	public Node(String name) throws IllegalArgumentException
	{
		this(name.toLowerCase(), name);
	}

	/**
	 * Create a new (non-virtual) node with a specific ID and name
	 * 
	 * @param ID the ID of this node - should be unique within a network!
	 * @param name the name of this node - will be used for displaying the node and is ideally unique, too
	 */
	public Node(String ID, String name) throws IllegalArgumentException
	{
		this(ID, name, false);
	}

	/**
	 * Create a new node with a specific ID, name and which is virtual or not
	 * 
	 * @param ID the ID of this node - should be unique within a network!
	 * @param name the name of this node - should be unique within a network!
	 * @param virtual whether this node is virtual
	 */
	public Node(String ID, String name, boolean virtual) throws IllegalArgumentException
	{
		if (ID == null)
		{
			String errormsg = "The ID of a node should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (name == null)
		{
			String errormsg = "The name of a node should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.ID = ID;
		this.name = name;
		this.virtual = virtual;
	}

	/**
	 * Return the unique ID of this node
	 * 
	 * @return the name of this node
	 */
	public String getID()
	{
		return ID;
	}

	/**
	 * Return the name of this node, in its original (unnormalized) form.
	 * 
	 * @return the name of this node
	 */
	public String getDisplayName()
	{
		return getDisplayName(false);
	}

	/**
	 * Return the name of this node, either normalized or in original form.
	 * 
	 * @param normalized whether or not the name should be converted to lower case.
	 * @return the name of this node
	 */
	public String getDisplayName(boolean normalized)
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
		String result = "Node " + ID;
		if (!ID.equals(name))
		{
			result += " (" + name + ")";
		}
		return result;
	}
}
