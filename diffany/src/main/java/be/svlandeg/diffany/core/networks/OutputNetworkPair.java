package be.svlandeg.diffany.core.networks;


/**
 * This class simply contains one differential network and its corresponding overlapping network.
 * 
 * @author Sofie Van Landeghem
 */
public class OutputNetworkPair
{
	
	private DifferentialNetwork dn;
	private ConsensusNetwork on;
	
	/**
	 * Create a new pair from one specific differential network and its corresponding overlap network. Neither of the two networks may be null!
     *
	 * @param dn the differential network
	 * @param on the overlap network
	 * @throws IllegalArgumentException when either of the provided networks is null
	 */
	public OutputNetworkPair(DifferentialNetwork dn, ConsensusNetwork on)
	{
		if (dn == null || on == null)
		{
			String errormsg = "The networks in the output network pair can not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.dn = dn;
		this.on = on;
	}
	
	/**
	 * Retrieve the differential network
	 * @return the differential network
	 */
	public DifferentialNetwork getDifferentialNetwork()
	{
		return dn;
	}
	
	/**
	 * Retrieve the overlapping network
	 * @return the overlapping network
	 */
	public ConsensusNetwork getOverlappingNetwork()
	{
		return on;
	}

}
