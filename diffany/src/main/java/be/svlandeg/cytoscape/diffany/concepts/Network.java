package be.svlandeg.cytoscape.diffany.concepts;

/**
 * Abstract class that represents a network: a collection of edges and nodes
 * @author Sofie Van Landeghem
 *
 */
public abstract class Network
{
	
	protected String name;  
	
	/**
	 * The name of this network should be enforced to be unique within one project.
	 * @param name the name of this network
	 */
	public Network(String name)
	{
		this.name = name;
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

}
