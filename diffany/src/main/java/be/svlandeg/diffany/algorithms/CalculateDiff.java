package be.svlandeg.diffany.algorithms;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.ReferenceNetwork;

/**
 * This class serves as an abstract layer between a GUI and the actual technical implementations 
 * of various algorithms that can calculate differential networks.
 * Depending on the parameters, this class will decide which algorithm to call.
 * 
 * @author Sofie Van Landeghem
 */
public class CalculateDiff 
{
	
	/**
	 * Calculate the differential network
	 * @param reference the reference network
	 * @param condition a condition-specific network
	 * @return the differential network between the two
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, ConditionNetwork condition)
	{
		throw new UnsupportedOperationException();
	}

}
