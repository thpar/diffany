package be.svlandeg.diffany.algorithms;

import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.concepts.*;
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
	
	protected boolean default_MIN = true;
	
	public double default_cutoff = 0.0;
	protected static String diffnameprefix = "diff_";
	protected static String overlapnameprefix = "overlap_";
	
	/**
	 * The constructor initializes the algorithm.
	 */
	public CalculateDiff()
	{}

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
	 * 	The overlapping network will get this name + the prefix 'overlap_'.
	 * @param cutoff the minimal value of a resulting edge for it to be included in the differential/overlapping networks
	 * @param log the logger that records logging messages
	 * 
	 * @return the differential network between the two
	 * @throws IllegalArgumentException if any of the parameters are null
	 */
	private DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, ConditionNetwork condition, EdgeOntology eo, 
			NodeMapper nm, String diff_name, double cutoff, Logger log) throws IllegalArgumentException
	{
		if (diff_name == null || reference == null || condition == null || eo == null || nm == null)
		{
			String errormsg = "Found null parameter in calculateDiffNetwork!";
			throw new IllegalArgumentException(errormsg);
		}
		Set<ConditionNetwork> allConditions = new HashSet<ConditionNetwork>();
		allConditions.add(condition);
		DifferentialNetwork dn = new CalculateDiffOfMore(log).calculateDiffNetwork(reference, allConditions, eo, nm, diff_name, cutoff);
		OverlappingNetwork sn = new CalculateDiffOfTwo(log).calculateOverlappingNetwork(reference, condition, eo, nm, overlapnameprefix + diff_name, cutoff, default_MIN);
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
	 * @param log the logger that records logging messages
	 * 
	 * @return the overlapping network between the two
	 * @throws IllegalArgumentException if any of the parameters are null
	 */
	public OverlappingNetwork calculateOverlappingNetwork(Network reference, Network condition, EdgeOntology eo, 
			NodeMapper nm, String overlap_name, double cutoff, Logger log) throws IllegalArgumentException
	{
		if (reference == null || condition == null || eo == null || nm == null || overlap_name == null)
		{
			String errormsg = "Found null parameter in calculateDiffNetwork!";
			throw new IllegalArgumentException(errormsg);
		}
		OverlappingNetwork sn = new CalculateDiffOfTwo(log).calculateOverlappingNetwork(reference, condition, eo, nm, overlap_name, cutoff, default_MIN);
		return sn;
	}

	/**
	 * Calculate all pairwise differential networks between the reference and each condition-specific network in the project.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object for each differential network.
	 * 
	 * The name of the differential network will be the name of the condition-specific network with prefix 'diff_'
	 * The name of the corresponding overlapping network will get an additional prefix 'overlap_'.
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * All calculated differential networks and the logger object are added to the project directly.
	 * The logger object will first be cleaned.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param configurationID the configuration ID of the configuration that needs to be run
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateAllPairwiseDifferentialNetworks(Project p, int configurationID) throws IllegalArgumentException
	{
		calculateAllPairwiseDifferentialNetworks(p, configurationID, default_cutoff);
	}
	
	/**
	 * Calculate all pairwise differential networks between the reference and each condition-specific network in the project.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object for each differential network.
	 * 
	 * The name of the differential network will be the name of the condition-specific network with prefix 'diff_'
	 * The name of the corresponding overlapping network will get an additional prefix 'overlap_'.
	 * 
	 * All calculated differential networks and the logger object are added to the project directly. 
	 * The logger object will first be cleaned.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param configurationID the configuration ID of the configuration that needs to be run
	 * @param cutoff the minimal value of a resulting edge for it to be included in the differential/overlapping networks
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateAllPairwiseDifferentialNetworks(Project p, int configurationID, double cutoff) throws IllegalArgumentException
	{
		EdgeOntology eo = p.getEdgeOntology();
		NodeMapper nm = p.getNodeMapper();
		Logger log = p.getLogger(configurationID);
		log.clean();
		
		RunConfiguration rc = p.getRunConfiguration(configurationID);
		ReferenceNetwork r = rc.getReferenceNetwork();
		for (ConditionNetwork c : rc.getConditionNetworks())
		{
			String diff_name = diffnameprefix + c.getName();
			log.log("Calculating the differential and overlap network between " + r.getName() + " and " + c.getName());
			DifferentialNetwork diff = calculateDiffNetwork(r, c, eo, nm, diff_name, cutoff, log);
			rc.addDifferential(diff);
		}
		log.log("Done!");
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
	 * @param log the logger that records logging messages
	 * 
	 * @return the overlapping network between all input networks
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public OverlappingNetwork calculateOverlappingNetwork(Set<Network> networks, EdgeOntology eo, 
			NodeMapper nm, String overlapping_name, double cutoff, Logger log) throws IllegalArgumentException
	{
		if (networks == null  || networks.isEmpty() || eo == null || nm == null || overlapping_name == null)
		{
			String errormsg = "Found null parameter in calculateOverlappingNetwork!";
			throw new IllegalArgumentException(errormsg);
		}
		if (networks.size() == 2)
		{
			OverlappingNetwork sn = new CalculateDiffOfTwo(log).calculateOverlappingNetwork(networks.iterator().next(), networks.iterator().next(), eo, 
					nm, overlapping_name, cutoff, default_MIN);
			return sn;
		}
		OverlappingNetwork sn = new CalculateDiffOfMore(log).calculateOverlappingNetwork(networks, eo, nm, overlapping_name, cutoff, default_MIN);
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
	 * 	The overlapping network will get this name + the prefix 'overlap_'.
	 * @param cutoff the minimal value of a resulting edge for it to be included in the differential/overlapping networks
	 * @param log the logger that records logging messages
	 * 
	 * @return the differential network between the two
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	private DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditions, EdgeOntology eo, 
			NodeMapper nm, String diff_name, double cutoff, Logger log) throws IllegalArgumentException
	{
		if (reference == null || conditions == null || conditions.isEmpty() || eo == null || nm == null || diff_name == null)
		{
			String errormsg = "The edge weight should be positive!";
			throw new IllegalArgumentException(errormsg);
		}
		DifferentialNetwork dn = new CalculateDiffOfMore(log).calculateDiffNetwork(reference, conditions, eo, nm, diff_name, cutoff);
		if (conditions.size() == 1)
		{
			OverlappingNetwork sn = new CalculateDiffOfTwo(log).calculateOverlappingNetwork(reference, conditions.iterator().next(), eo, 
					nm, overlapnameprefix + diff_name, cutoff, default_MIN);
			dn.setOverlappingNetwork(sn);
			return dn;
		}
		Set<Network> allNetworks = new HashSet<Network>();
		allNetworks.add(reference);
		allNetworks.addAll(conditions);
		OverlappingNetwork sn = new CalculateDiffOfMore(log).calculateOverlappingNetwork(allNetworks, eo, nm, overlapnameprefix + diff_name, cutoff, default_MIN);
		dn.setOverlappingNetwork(sn);
		return dn;
	}
	

	/**
	 * Calculate the differential network between the reference and a set of condition-specific networks.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
	 * 
	 * The name of the differential network will be 'all_conditions_against_reference_diff'.
	 * The name of the corresponding overlapping network will get an additional prefix 'overlap_' at the end.
	 * 
	 * The calculated differential networks and the logger object are added are added to the project directly.
	 * The logger object will first be cleaned.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param configurationID the configuration ID of the configuration that needs to be run
	 * @param diff_name the name to give to the differential network. 
	 * 	The overlapping network will get this name + the prefix 'overlap_'.
	 * @param cutoff the minimal value of a resulting edge for it to be included in the overlapping network
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateOneDifferentialNetwork(Project p, int configurationID, String diff_name, double cutoff) throws IllegalArgumentException
	{
		EdgeOntology eo = p.getEdgeOntology();
		NodeMapper nm = p.getNodeMapper();
		Logger log = p.getLogger(configurationID);
		log.clean();
		
		RunConfiguration rc = p.getRunConfiguration(configurationID);
		ReferenceNetwork r = rc.getReferenceNetwork();
		
		Set<ConditionNetwork> cs = new HashSet<ConditionNetwork>(rc.getConditionNetworks());
		log.log("Calculating the differential and overlap network between " + r.getName() + " and " 
				+ cs.size() + " condition-dependent network(s)");
		DifferentialNetwork diff = calculateDiffNetwork(r, cs, eo, nm, diff_name, cutoff, log);
		rc.addDifferential(diff);
		log.log("Done!");
	}
	
	/**
	 * Calculate the differential network between the reference and a set of condition-specific networks.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
	 * 
	 * The name of the differential network will be 'all_conditions_against_reference_diff'.
	 * The name of the corresponding overlapping network will get an additional prefix 'overlap_' at the end.
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * The calculated differential networks and the logger object are added are added to the project directly.
	 * The logger object will first be cleaned.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param configurationID the configuration ID of the configuration that needs to be ru
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateOneDifferentialNetwork(Project p, int configurationID) throws IllegalArgumentException
	{
		String diff_name = diffnameprefix + "all_conditions_against_reference";
		calculateOneDifferentialNetwork(p, configurationID, diff_name, default_cutoff);
	}
	
	/**
	 * Calculate the differential network between the reference and a set of condition-specific networks.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
	 * 
	 * The name of the differential network will be 'all_conditions_against_reference_diff'.
	 * The name of the corresponding overlapping network will get an additional prefix 'overlap_' at the end.
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * The calculated differential networks and the logger object are added are added to the project directly.
	 * The logger object will first be cleaned.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param cutoff the minimal value of a resulting edge for it to be included in the overlapping network
	 * @param configurationID the configuration ID of the configuration that needs to be ru
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateOneDifferentialNetwork(Project p, int configurationID, double cutoff) throws IllegalArgumentException
	{
		String diff_name = diffnameprefix + "all_conditions_against_reference";
		calculateOneDifferentialNetwork(p, configurationID, diff_name, cutoff);
	}
	
	/**
	 * Calculate the differential network between the reference and a set of condition-specific networks.
	 * Adds the 'overlapping' or 'house-keeping' interactions as a related OverlappingNetwork object.
	 * 
	 * The name of the differential network will be 'all_conditions_against_reference_diff'.
	 * The name of the corresponding overlapping network will get an additional prefix 'overlap_' at the end.
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * The calculated differential networks and the logger object are added are added to the project directly.
	 * The logger object will first be cleaned.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param configurationID the configuration ID of the configuration that needs to be run
	 * @param diff_name the name to give to the differential network. 
	 * 	The overlapping network will get this name + the prefix 'overlap_'.
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateOneDifferentialNetwork(Project p, int configurationID, String diff_name) throws IllegalArgumentException
	{
		calculateOneDifferentialNetwork(p, configurationID, diff_name, default_cutoff);
	}
}
