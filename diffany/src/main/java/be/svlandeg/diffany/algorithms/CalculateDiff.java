package be.svlandeg.diffany.algorithms;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * This class serves as an abstract layer between a GUI and the actual technical
 * implementations of various algorithms that can calculate differential
 * networks. Depending on the parameters, this class will decide which algorithm
 * to call.
 * 
 * @author Sofie Van Landeghem
 */
public class CalculateDiff
{

	/**
	 * Calculate the differential network between the reference and
	 * condition-specific network.
	 * 
	 * @param reference
	 *            the reference network
	 * @param condition
	 *            a condition-specific network
	 * @param eo
	 *            the edge ontology that provides meaning to the edge types
	 * @param nm
	 *            the node mapper that allows to map nodes from the one network
	 *            to the other
	 * @param diff_name
	 *            the name to give to the differential network
	 * @return the differential network between the two
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, ConditionNetwork condition, EdgeOntology eo, NodeMapper nm, String diff_name)
	{
		return new CalculateDiffOfTwo().calculateDiffNetwork(reference, condition, eo, nm, diff_name);
	}

}
