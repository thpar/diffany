package be.svlandeg.diffany.algorithms;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * This class serves as an abstract layer between a GUI and the actual technical implementations of various algorithms 
 * that can calculate differential networks. 
 * Depending on the parameters, this class will decide which algorithm to call.
 * 
 * @author Sofie Van Landeghem
 */
public class CalculateDiff
{
	
	protected static String diffnamesuffix = "_diff";
	protected static String overlapnamesuffix = "_overlap";

	/**
	 * Calculate the differential network between the reference and
	 * condition-specific network.
	 * 
	 * @param reference the reference network
	 * @param condition a condition-specific network
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network. 
	 * 	The shared network will get this name + the appendix '_overlap'.
	 * @return the differential network between the two
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, ConditionNetwork condition, EdgeOntology eo, NodeMapper nm, String diff_name)
	{
		return new CalculateDiffOfTwo().calculateDiffNetwork(reference, condition, eo, nm, diff_name);
	}
	
	/**
	 * Calculate the differential network between the reference and condition-specific network. 
	 * The name of the differential network will be the name of the condition-specific network with appendix '_diff'.
	 * The name of the corresponding shared network will get an additional appendix '_overlap'.
	 * 
	 * @param reference the reference network
	 * @param condition a condition-specific network
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @return the differential network between the two
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, ConditionNetwork condition, EdgeOntology eo, NodeMapper nm)
	{
		String diff_name = condition.getName() + diffnamesuffix;
		return new CalculateDiffOfTwo().calculateDiffNetwork(reference, condition, eo, nm, diff_name);
	}

	/**
	 * Calculate all pairwise differential networks between the reference and each condition-specific network in the project.
	 * The name of the differential network will be the name of the condition-specific network with appendix '_diff'
	 * The name of the corresponding shared network will get an additional appendix '_overlap'.
	 * 
	 * All calculated differential networks are added to the project directly.
	 * 
	 * @param p the project which stores the reference and condition-specific networks.
	 */
	public void calculateAllPairwiseDifferentialNetworks(Project p)
	{
		ReferenceNetwork r = p.getReferenceNetwork();
		EdgeOntology eo = p.getEdgeOntology();
		NodeMapper nm = p.getNodeMapper();
		for (ConditionNetwork c : p.getConditionNetworks())
		{
			DifferentialNetwork diff = calculateDiffNetwork(r, c, eo, nm);
			p.addDifferential(diff);
		}
	}

}
