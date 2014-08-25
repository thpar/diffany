package be.svlandeg.diffany.core.networks;


/**
 * This class simply contains one differential network and its corresponding overlapping network.
 * 
 * @author Sofie Van Landeghem
 */
public class OutputNetworkPair
{
	
	private DifferentialNetwork dn;
	private ConsensusNetwork cn;
	
	/**
	 * Create a new pair from one specific differential network and its counterpart consensus network. Neither of the two networks may be null!
     *
	 * @param dn the differential network
	 * @param cn the consensus network
	 * @throws IllegalArgumentException when either of the provided networks is null
	 */
	public OutputNetworkPair(DifferentialNetwork dn, ConsensusNetwork cn)
	{
		if (dn == null || cn == null)
		{
			String errormsg = "The networks in the output network pair can not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.dn = dn;
		this.cn = cn;
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
	 * Retrieve the consensus network
	 * @return the consensus network
	 */
	public ConsensusNetwork getConsensusNetwork()
	{
		return cn;
	}

}
