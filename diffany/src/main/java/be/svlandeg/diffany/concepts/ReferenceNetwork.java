package be.svlandeg.diffany.concepts;


/**
 * A reference network is used as comparison against condition-dependent networks.
 * It can be seen as a 'static' network with unspecified/unknown/wild-type conditions.
 * 
 * @author Sofie Van Landeghem
 *
 */
public class ReferenceNetwork extends Network
{
	
	protected String name;
	
	/**
	 * Create a new static reference network.
	 * 
	 * @param name the name of this network
	 * 
	 */
	public ReferenceNetwork(String name)
	{
		super(name);
	}

	/* (non-Javadoc)
	 * @see be.svlandeg.diffany.concepts.Network#getStringRepresentation()
	 */
	@Override
	public String getStringRepresentation()
	{
		return name + ": static reference network";
	}

	

	
}
