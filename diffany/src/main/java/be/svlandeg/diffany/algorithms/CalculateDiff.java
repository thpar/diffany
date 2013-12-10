package be.svlandeg.diffany.algorithms;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.concepts.OverlappingNetwork;
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
	
	protected double default_cutoff = 0.0;
	protected static String diffnamesuffix = "_diff";
	protected static String overlapnamesuffix = "_overlap";
	
	protected CalculateDiffOfTwo twoProcessor;
	protected CalculateDiffOfMore moreProcessor;
	
	public CalculateDiff()
	{
		twoProcessor = new CalculateDiffOfTwo();
		moreProcessor = new CalculateDiffOfMore();
	}

	/////// METHODS BETWEEN 1 REFERENCE NETWORK AND 1 CONDITION-SPECIFIC NETWORK  //////////////
	
	/**
	 * Calculate the differential network between the reference and condition-specific network. 
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
	 * 
	 * @param reference the reference network
	 * @param condition a condition-specific network
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network. 
	 * 	The overlapping network will get this name + the appendix '_overlap'.
	 * @param cutoff the minimal value of a resulting edge for it to be included in the differential/overlapping networks
	 * 
	 * @return the differential network between the two
	 * @throws IllegalArgumentException if any of the parameters are null
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, ConditionNetwork condition, EdgeOntology eo, 
			NodeMapper nm, String diff_name, double cutoff) throws IllegalArgumentException
	{
		if (reference == null || condition == null || eo == null || nm == null || diff_name == null)
		{
			String errormsg = "Found null parameter in calculateDiffNetwork!";
			throw new IllegalArgumentException(errormsg);
		}
		DifferentialNetwork dn =  twoProcessor.calculateDiffNetwork(reference, condition, eo, nm, diff_name, cutoff);
		OverlappingNetwork sn = twoProcessor.calculateOverlappingNetwork(reference, condition, eo, nm, diff_name + overlapnamesuffix, cutoff);
		dn.setOverlappingNetwork(sn);
		return dn;
	}
	
	/**
	 * Calculate the overlapping network between the reference and condition-specific network, which includes the 'house-keeping' interactions 
	 * 
	 * @param reference the reference network
	 * @param condition a condition-specific network
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param overlap_name the name to give to the overlapping network
	 * @param cutoff the minimal value of a resulting edge for it to be included in the differential/overlapping networks
	 * 
	 * @return the overlapping network between the two
	 * @throws IllegalArgumentException if any of the parameters are null
	 */
	public OverlappingNetwork calculateOverlappingNetwork(Network reference, Network condition, EdgeOntology eo, 
			NodeMapper nm, String overlap_name, double cutoff) throws IllegalArgumentException
	{
		if (reference == null || condition == null || eo == null || nm == null || overlap_name == null)
		{
			String errormsg = "Found null parameter in calculateDiffNetwork!";
			throw new IllegalArgumentException(errormsg);
		}
		OverlappingNetwork sn = twoProcessor.calculateOverlappingNetwork(reference, condition, eo, nm, overlap_name, cutoff);
		return sn;
	}
	
	/**
	 * Calculate the differential network between the reference and condition-specific network.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
	 * 
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * @param reference the reference network
	 * @param condition a condition-specific network
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network. 
	 * 	The overlapping network will get this name + the appendix '_overlap'.
	 * 
	 * @return the differential network between the two
	 * @throws IllegalArgumentException if any of the parameters are null
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, ConditionNetwork condition, EdgeOntology eo, 
			NodeMapper nm, String diff_name) throws IllegalArgumentException
	{
		return calculateDiffNetwork(reference, condition, eo, nm, diff_name, default_cutoff);
	}
	
	/**
	 * Calculate the differential network between the reference and condition-specific network. 
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
	 * 
	 * The name of the differential network will be the name of the condition-specific network with appendix '_diff'.
	 * The name of the corresponding overlapping network will get an additional appendix '_overlap'.
	 * 
	 * @param reference the reference network
	 * @param condition a condition-specific network
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param cutoff the minimal value of a resulting edge for it to be included in the differential/overlapping networks
	 * 
	 * @return the differential network between the two
	 * @throws IllegalArgumentException if any of the parameters are null
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, ConditionNetwork condition, EdgeOntology eo, 
			NodeMapper nm, double cutoff) throws IllegalArgumentException
	{
		String diff_name = condition.getName() + diffnamesuffix;
		return calculateDiffNetwork(reference, condition, eo, nm, diff_name, cutoff);
	}
	
	/**
	 * Calculate the differential network between the reference and condition-specific network. 
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
	 * 
	 * The name of the differential network will be the name of the condition-specific network with appendix '_diff'.
	 * The name of the corresponding overlapping network will get an additional appendix '_overlap'.
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * @param reference the reference network
	 * @param condition a condition-specific network
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * 
	 * @return the differential network between the two
	 * @throws IllegalArgumentException if any of the parameters are null
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, ConditionNetwork condition, EdgeOntology eo, 
			NodeMapper nm) throws IllegalArgumentException
	{
		String diff_name = condition.getName() + diffnamesuffix;
		return calculateDiffNetwork(reference, condition, eo, nm, diff_name);
	}

	/**
	 * Calculate all pairwise differential networks between the reference and each condition-specific network in the project.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object for each differential network.
	 * 
	 * The name of the differential network will be the name of the condition-specific network with appendix '_diff'
	 * The name of the corresponding overlapping network will get an additional appendix '_overlap'.
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * All calculated differential networks are added to the project directly.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateAllPairwiseDifferentialNetworks(Project p) throws IllegalArgumentException
	{
		calculateAllPairwiseDifferentialNetworks(p, default_cutoff);
	}
	
	/**
	 * Calculate all pairwise differential networks between the reference and each condition-specific network in the project.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object for each differential network.
	 * 
	 * The name of the differential network will be the name of the condition-specific network with appendix '_diff'
	 * The name of the corresponding overlapping network will get an additional appendix '_overlap'.
	 * 
	 * All calculated differential networks are added to the project directly.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param cutoff the minimal value of a resulting edge for it to be included in the differential/overlapping networks
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateAllPairwiseDifferentialNetworks(Project p, double cutoff) throws IllegalArgumentException
	{
		ReferenceNetwork r = p.getReferenceNetwork();
		EdgeOntology eo = p.getEdgeOntology();
		NodeMapper nm = p.getNodeMapper();
		for (ConditionNetwork c : p.getConditionNetworks())
		{
			DifferentialNetwork diff = calculateDiffNetwork(r, c, eo, nm, cutoff);
			p.addDifferential(diff);
		}
	}
	
	
	/////// METHODS BETWEEN 1 REFERENCE NETWORK AND MULTIPLE CONDITION-SPECIFIC NETWORKS  //////////////
	
	/**
	 * Calculate the overlapping network for a set of input networks (at least 2).
	 * 
	 * @param networks a set of networks (at least 1)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param overlapping_name the name to give to the overlapping network. 
	 * @param cutoff the minimal value of a resulting edge for it to be included in the overlapping network
	 * 
	 * @return the overlapping network between all input networks
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public OverlappingNetwork calculateOverlappingNetwork(Set<Network> networks, EdgeOntology eo, 
			NodeMapper nm, String overlapping_name, double cutoff) throws IllegalArgumentException
	{
		if (networks == null  || networks.isEmpty() || eo == null || nm == null || overlapping_name == null)
		{
			String errormsg = "Found null parameter in calculateOverlappingNetwork!";
			throw new IllegalArgumentException(errormsg);
		}
		if (networks.size() == 2)
		{
			OverlappingNetwork sn = twoProcessor.calculateOverlappingNetwork(networks.iterator().next(), networks.iterator().next(), eo, nm, overlapping_name, cutoff);
			return sn;
		}
		OverlappingNetwork sn = moreProcessor.calculateOverlappingNetwork(networks, eo, nm, overlapping_name, cutoff);
		return sn;
	}
	
	/**
	 * Calculate the differential network between the reference and a set of condition-specific networks.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
	 * 
	 * @param reference the reference network
	 * @param conditions a set of condition-specific networks (at least 1)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network. 
	 * 	The overlapping network will get this name + the appendix '_overlap'.
	 * @param cutoff the minimal value of a resulting edge for it to be included in the differential/overlapping networks
	 * 
	 * @return the differential network between the two
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditions, EdgeOntology eo, 
			NodeMapper nm, String diff_name, double cutoff) throws IllegalArgumentException
	{
		if (reference == null || conditions == null || conditions.isEmpty() || eo == null || nm == null || diff_name == null)
		{
			String errormsg = "The edge weight should be positive!";
			throw new IllegalArgumentException(errormsg);
		}
		if (conditions.size() == 1)
		{
			DifferentialNetwork dn = twoProcessor.calculateDiffNetwork(reference, conditions.iterator().next(), eo, nm, diff_name, cutoff);
			OverlappingNetwork sn = twoProcessor.calculateOverlappingNetwork(reference, conditions.iterator().next(), eo, nm, diff_name + overlapnamesuffix, cutoff);
			dn.setOverlappingNetwork(sn);
			return dn;
		}
		DifferentialNetwork dn = moreProcessor.calculateDiffNetwork(reference, conditions, eo, nm, diff_name, cutoff);
		Set<Network> allNetworks = new HashSet<Network>();
		allNetworks.add(reference);
		allNetworks.addAll(conditions);
		OverlappingNetwork sn = moreProcessor.calculateOverlappingNetwork(allNetworks, eo, nm, diff_name + overlapnamesuffix, cutoff);
		dn.setOverlappingNetwork(sn);
		return dn;
	}
	
	/**
	 * Calculate the differential network between the reference and a set of condition-specific networks.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
	 * 
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * @param reference the reference network
	 * @param conditions a set of condition-specific networks (at least 1)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network. 
	 * 	The overlapping network will get this name + the appendix '_overlap'.
	 * 
	 * @return the differential network between the two
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditions, EdgeOntology eo, 
			NodeMapper nm, String diff_name) throws IllegalArgumentException
	{
		return calculateDiffNetwork(reference, conditions, eo, nm, diff_name, default_cutoff);
	}
	
	/**
	 * Calculate the differential network between the reference and a set of condition-specific networks.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
	 * 
	 * The name of the differential network will be 'all_conditions_against_reference_diff'.
	 * The name of the corresponding overlapping network will be the same, but with '_overlap' at the end.
	 * 
	 * @param reference the reference network
	 * @param conditions a set of condition-specific networks (at least 1)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param cutoff the minimal value of a resulting edge for it to be included in the differential/overlapping networks
	 * 
	 * @return the differential network between the two
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditions, EdgeOntology eo, 
			NodeMapper nm, double cutoff) throws IllegalArgumentException
	{
		String diff_name = "all_conditions_against_reference" + diffnamesuffix;
		return calculateDiffNetwork(reference, conditions, eo, nm, diff_name, cutoff);
	}
	
	/**
	 * Calculate the differential network between the reference and a set of condition-specific networks.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
	 * 
	 * The name of the differential network will be 'all_conditions_against_reference_diff'.
	 * The name of the corresponding overlapping network will be the same, but with '_overlap' at the end.
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * @param reference the reference network
	 * @param conditions a set of condition-specific networks (at least 1)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * 
	 * @return the differential network between the two
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditions, EdgeOntology eo, 
			NodeMapper nm) throws IllegalArgumentException
	{
		String diff_name = "all_conditions_against_reference" + diffnamesuffix;
		return calculateDiffNetwork(reference, conditions, eo, nm, diff_name);
	}
	
	/**
	 * Calculate the differential network between the reference and all condition-specific networks in the project.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object for each differential network.
	 * 
	 * The name of the differential network will be 'all_conditions_against_reference_diff'.
	 * The name of the corresponding overlapping network will get an additional appendix '_overlap' at the end.
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * The calculated differential networks is added to the project directly.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateOneDifferentialNetwork(Project p) throws IllegalArgumentException
	{
		calculateOneDifferentialNetwork(p, default_cutoff);
	}
	
	/**
	 * Calculate all pairwise differential networks between the reference and each condition-specific network in the project.
	 * The name of the differential network will be the name of the condition-specific network with appendix '_diff'
	 * The name of the corresponding overlapping network will get an additional appendix '_overlap'.
	 * 
	 * The calculated differential networks is added to the project directly.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param cutoff the minimal value of a resulting edge for it to be included in the overlapping network
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateOneDifferentialNetwork(Project p, double cutoff) throws IllegalArgumentException
	{
		ReferenceNetwork r = p.getReferenceNetwork();
		EdgeOntology eo = p.getEdgeOntology();
		NodeMapper nm = p.getNodeMapper();
		Set<ConditionNetwork> cs = new HashSet<ConditionNetwork>(p.getConditionNetworks());
		DifferentialNetwork diff = calculateDiffNetwork(r, cs, eo, nm);
		p.addDifferential(diff);
	}

}
