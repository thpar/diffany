package be.svlandeg.diffany.algorithms;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import be.svlandeg.diffany.concepts.Condition;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.concepts.OverlappingNetwork;
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
	
	protected CalculateDiffOfTwo twoProcessor;
	
	/**
	 * Constructor initializes the algorithm suites.
	 */
	public CalculateDiffOfMore()
	{
		twoProcessor = new CalculateDiffOfTwo();
	}
	
	/**
	 * Calculate the differential network between the reference and condition-specific networks. 
	 * The overlapping network should be calculated independently!
	 * This method can only be called from within the package (CalculateDiff) and can thus assume proper input.
	 * 
	 * @param reference the reference network
	 * @param conditions a set of condition-specific networks (at least 2)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network. 
	 * 
	 * @return the differential network between the two    
	 */
	protected DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditionNetworks, EdgeOntology eo, 
			NodeMapper nm, String diff_name, double cutoff)
	{
		Set<Network> networks = new HashSet<Network>();
		networks.addAll(conditionNetworks);
		OverlappingNetwork onMIN = calculateOverlappingNetwork(networks, eo, nm, diff_name + "_overlap", cutoff, true);
		OverlappingNetwork onMAX = calculateOverlappingNetwork(networks, eo, nm, diff_name + "_overlap", cutoff, false);

		// TODO
		return null;
	}
	
	
	/**
	 * Calculate the overlapping network between a set of networks. 
	 * This method can only be called from within the package (CalculateDiff) and can thus assume proper input.
	 * 
	 * @param networks a set of networks (at least 2)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param overlapping_name the name to give to the overlapping network. 
	 * @param minOperator whether or not to take the minimum of the edge weights for the overlapping edges - if false, the maximum is taken
	 * 
	 * @return the differential network between the two
	 */
	protected OverlappingNetwork calculateOverlappingNetwork(Set<Network> networks, EdgeOntology eo, 
			NodeMapper nm, String overlapping_name, double cutoff, boolean minOperator)
	{
		List<Network> listedNetworks = new ArrayList<Network>();
		listedNetworks.addAll(networks);
		
		int numberOfNetworks = listedNetworks.size();
		int first = 0;
		int second = 1;
		
		
		Network firstN = listedNetworks.get(first);
		Network secondN = listedNetworks.get(second);
		OverlappingNetwork overlapTmp = twoProcessor.calculateOverlappingNetwork(firstN, secondN, eo, nm, overlapping_name, cutoff, minOperator);
		second++;
		
		while (second < numberOfNetworks)
		{	
			secondN = listedNetworks.get(second);
			overlapTmp = twoProcessor.calculateOverlappingNetwork(overlapTmp, secondN, eo, nm, overlapping_name, cutoff, minOperator);
			second++;
		}
		
		return overlapTmp;
	}
	
}
