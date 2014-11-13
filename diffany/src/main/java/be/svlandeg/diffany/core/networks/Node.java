package be.svlandeg.diffany.core.networks;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;


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
	protected Map<String, String> attributes;
	
	// TODO - record this information somewhere else?
	public static final String de_attribute = "differentially_expressed";
	public static final String phos_attribute = "phosphorylation_site";
	public static final String kinase_attribute = "kinase_function";
	
	public static final String upregulated = "up-regulated";
	public static final String downregulated = "down-regulated";
	public static final String not_de = "no";

	/**
	 * Create a new node with a specific ID, name and which is virtual or not
	 * 
	 * @param ID the ID of this node - should be unique within a network!
	 * @param name the name of this node - should be unique within a network!
	 * @param virtual whether this node is virtual
	 * @throws IllegalArgumentException when the name or ID is null
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
		attributes = new HashMap<String, String>();
	}
	
	/**
	 * Create a new (non-virtual) node with a specific ID and name
	 * 
	 * @param ID the ID of this node - should be unique within a network!
	 * @param name the name of this node - will be used for displaying the node and is ideally unique, too
	 * @throws IllegalArgumentException when the name or ID is null
	 */
	public Node(String ID, String name) throws IllegalArgumentException
	{
		this(ID, name, false);
	}
	
	/**
	 * Create a new (non-virtual) node with a specific name. The lower-case version of this name will be used as unique identifier, so ensure its unambiguity across the project!
	 * The use of this constructor should be avoided!
	 * 
	 * @param name the name of this node - should be unique within a network!
	 * @throws IllegalArgumentException when the name is null
	 */
	public Node(String name) throws IllegalArgumentException
	{
		this(name.toLowerCase(), name);
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
	 * Retrieve all attributes recorded for this node.
	 * Note that from the Network this Node belongs to, you can derive all required attribute names for this node.
	 * 
	 * @return all attribute names recorded for this Node
	 */
	public Set<String> getAllAttributeNames()
	{
		return attributes.keySet();
	}
	
	
	/**
	 * Retrieve the value of a certain attribute, or null when it is not defined for this node.
	 * From the Network this Node belongs to, you can derive all required attribute names for this node.
	 * 
	 * @param attributeName the name of the attribute
	 * @return the value of the attribute
	 */
	public String getAttribute(String attributeName)
	{
		return attributes.get(attributeName);
	}
	
	/**
	 * Set the value of a certain attribute. If a value already existed for this attribute, it is overwritten.
	 * 
	 * @param attributeName the name of the attribute (should not be null!)
	 * @param value the value of the attribute (should not be null!)
	 * 
	 * @throws IllegalArgumentException when either of the two input parameters is null
	 */
	public void setAttribute(String attributeName, String value)
	{
		if (attributeName == null)
		{
			String errormsg = "The attribute name should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (name == null)
		{
			String errormsg = "The attribute value should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		attributes.put(attributeName, value);
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
	
	/**
	 * Provide a long string of this node, including its attributes
	 * @return a full string representation of this node
	 */
	public String toLongString()
	{
		String result = "Node " + ID;
		if (!ID.equals(name))
		{
			result += " (" + name + ")";
		}
		for (String att : attributes.keySet())
		{
			result += '\t' + att + "=" + attributes.get(att);
		}
		return result;
	}
}
