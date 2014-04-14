package be.svlandeg.diffany.core.project;

import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.OverlappingNetwork;

/**
 * This class keeps the output of differential network algorithms, both differential networks as well overlapping networks.
 * Either of the two can be null when it was not calculated.
 * 
 * @author Sofie Van Landeghem
 */
public class DifferentialOutput
{
	
	private DifferentialNetwork dn;
	private OverlappingNetwork on;

	
	/**
	 * Create a new empty output object from both a non-null differential and a non-null overlapping network.
     *
	 * @param dn the differential network
	 * @param on the overlap network
	 * @throws IllegalArgumentException when either of the provided networks is null
	 */
	public DifferentialOutput(DifferentialNetwork dn, OverlappingNetwork on)
	{
		setDifferential(dn);
		setOverlap(on);
	}
	
	/**
	 * Create a new empty output object.
     *
	 * @throws IllegalArgumentException when either of the provided networks is null
	 */
	public DifferentialOutput()
	{
		this.dn = null;
		this.on = null;
	}
	
	/**
	 * Set the differential output network.
	 * @param dn the differential output network
	 */
	public void setDifferential(DifferentialNetwork dn)
	{
		if (dn == null)
		{
			String errormsg = "The specified output differential network(s) can not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.dn = dn;
	}
	
	/**
	 * Set the overlap output network.
	 * @param on the overlap output network
	 */
	public void setOverlap(OverlappingNetwork on)
	{
		if (on == null)
		{
			String errormsg = "The specified output overlap network(s) can not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.on = on;
	}
	
	/**
	 * Retrieve the differential and overlap networks as a result pair
	 * @return the result pair, containing both the differential and overlap networks
	 * @throws IllegalArgumentException when either of the two was null
	 */
	public OutputNetworkPair getOutputAsPair()
	{
		if (dn == null)
		{
			String errormsg = "Can not provided the differential network: it was not calculated!";
			throw new IllegalArgumentException(errormsg);
		}
		if (on == null)
		{
			String errormsg = "Can not provided the overlapping network: it was not calculated!";
			throw new IllegalArgumentException(errormsg);
		}
		return new OutputNetworkPair(dn, on);
	}
	
	/**
	 * Retrieve the differential network
	 * @return the differential network
	 * @throws IllegalArgumentException when the differential network was null (i.e. not calculated)
	 */
	public DifferentialNetwork getDifferentialNetwork()
	{
		if (dn == null)
		{
			String errormsg = "Can not provided the differential network: it was not calculated!";
			throw new IllegalArgumentException(errormsg);
		}
		return dn;
	}
	
	/**
	 * Retrieve the overlapping network
	 * @return the overlapping network
	 * @throws IllegalArgumentException when the overlapping network was null (i.e. not calculated)
	 */
	public OverlappingNetwork getOverlappingNetwork()
	{
		if (on == null)
		{
			String errormsg = "Can not provided the overlapping network: it was not calculated!";
			throw new IllegalArgumentException(errormsg);
		}
		return on;
	}
	
}
