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
	protected Map<String, Object> attributes;
	
	// TODO v2.2: record this information somewhere else
	public static final String de_attribute = "differentially_expressed";
	public static final String phos_attribute = "phosphorylation_site";
	public static final String kinase_attribute = "kinase_function";
	
	public static final String upregulated = "up-regulated";
	public static final String downregulated = "down-regulated";
	public static final String not_de = "no";

	
	/**
	 * Create a new (non-virtual) node with a specific ID and name
	 * 
	 * @param ID the ID of this node - should be unique within a network!
	 * @param name the name of this node - will be used for displaying the node and is ideally unique, too
	 * @throws IllegalArgumentException when the name or ID is null
	 */
	public Node(String ID, String name) throws IllegalArgumentException
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
		attributes = new HashMap<String, Object>();
	}

	/**
	 * Return the unique ID of this node.
	 * 
	 * @return the name of this node
	 */
	public String getID()
	{
		return ID;
	}

	/**
	 * Return the name of this node.
	 * 
	 * @return the name of this node
	 */
	public String getDisplayName()
	{
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
	public Object getAttribute(String attributeName)
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
	public void setAttribute(String attributeName, Object value)
	{
		if (attributeName == null)
		{
			String errormsg = "The attribute name should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (name == null)
		{
			//FIXME: is this checking the name and not the attribute value on purpose?
			//Attributes SHOULD be allowed to be null though...
			String errormsg = "The attribute value should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		attributes.put(attributeName, value);
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
