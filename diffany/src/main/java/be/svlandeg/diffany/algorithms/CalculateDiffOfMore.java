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
	
	
	/**
	 * Calculate the differential network between the reference and condition-specific networks. 
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
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
		OverlappingNetwork sn = calculateOverlappingNetwork(networks, eo, nm, diff_name + "_overlap", cutoff);
		
		Set<Condition> conditions = new HashSet<Condition>();
		for (ConditionNetwork c : conditionNetworks)
		{
			conditions.addAll(c.getConditions());
		}
		ConditionNetwork cn = new ConditionNetwork("All_conditions", conditions);
		cn.setNodesAndEdges(sn.getNodes(), sn.getEdges());
		
		DifferentialNetwork dn = new CalculateDiff().calculateDiffNetwork(reference, cn, eo, nm, diff_name, cutoff);
		dn.setOverlappingNetwork(sn);
		
		return dn;
	}
	
	/**
	 * Calculate the overlapping network between a set of networks. 
	 * This method can only be called from within the package (CalculateDiff) and can thus assume proper input.
	 * 
	 * @param networks a set of networks (at least 2)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param overlapping_name the name to give to the overlapping network. 
	 * 
	 * @return the differential network between the two
	 */
	protected OverlappingNetwork calculateOverlappingNetwork(Set<Network> networks, EdgeOntology eo, 
			NodeMapper nm, String overlapping_name, double cutoff)
	{
		List<Network> listedNetworks = new ArrayList<Network>();
		listedNetworks.addAll(networks);
		
		int numberOfNetworks = listedNetworks.size();
		int first = 0;
		int second = 1;
		
		
		Network firstN = listedNetworks.get(first);
		Network secondN = listedNetworks.get(second);
		OverlappingNetwork overlapTmp = new CalculateDiff().calculateOverlappingNetwork(firstN, secondN, eo, nm, overlapping_name, cutoff);
		second++;
		
		while (second < numberOfNetworks)
		{	
			secondN = listedNetworks.get(second);
			overlapTmp = new CalculateDiff().calculateOverlappingNetwork(overlapTmp, secondN, eo, nm, overlapping_name, cutoff);
			second++;
		}
		
		return overlapTmp;
	}
	
}
