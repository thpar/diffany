package be.svlandeg.diffany.algorithms;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * This class can calculate a differential network between one reference and a set of condition-specific networks. 
 * Currently this algorithm assumes a 1-to-1 mapping of nodes between the two networks!
 * 
 * @author Sofie Van Landeghem
 */
public class CalculateDiffOfMore
{
	
	
	/**
	 * Calculate the differential network between the reference and condition-specific networks. 
	 * The differential network will also store the shared/overlap/'house-keeping' interactions.
	 * This method can only be called from within the package (CalculateDiff) and can thus assume proper input.
	 * 
	 * @param reference the reference network
	 * @param conditions a set of condition-specific networks (at least 2)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network. 
	 * 	The associated shared network will be named similarly with appendix '_overlap'.
	 * @return the differential network between the two
	 *         
	 * TODO: expand this algorithm to be able to deal with n-m node mappings (v.2.0) 
	 * TODO: expand this algorithm to be able to deal with more than 1 edge between two nodes 
	 * in the original networks (v.1.0)
	 */
	protected DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditions, EdgeOntology eo, 
			NodeMapper nm, String diff_name, double cutoff)
	{
		List<ConditionNetwork> listedConditions = new ArrayList<ConditionNetwork>();
		listedConditions.addAll(conditions);
		
		// TODO
		
		return null;
	}
	

}
