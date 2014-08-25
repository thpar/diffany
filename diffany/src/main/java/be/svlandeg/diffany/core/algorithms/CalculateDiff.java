package be.svlandeg.diffany.core.algorithms;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunConfiguration;
import be.svlandeg.diffany.core.project.RunDiffConfiguration;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

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

	public double default_weight_cutoff = 0.0;
	protected static String diffnameprefix = "diff_";
	protected static String consensusnameprefix = "consensus_";

	public enum RunMode
	{
		EDGEBYEDGE, MACHINELEARNING
	};

	protected RunMode mode;

	/**
	 * Initialize the algorithm. By default, edge-by-edge comparisons are selected
	 */
	public CalculateDiff()
	{
		mode = RunMode.EDGEBYEDGE;
	}

	/**
	 * Initialize the algorithm, depending on the mode parameter.
	 * 
	 * @param mode the mode with which consensus/differential algorithms will be calculated.
	 */
	public CalculateDiff(RunMode mode)
	{
		this.mode = mode;
	}

	/////// METHODS BETWEEN 1 REFERENCE NETWORK AND MULTIPLE CONDITION-SPECIFIC NETWORKS  //////////////

	/**
	 * Calculate the consensus network for a set of input networks (at least 2).
	 * 
	 * @param networks a set of networks (at least 1)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param consensus_name the name to give to the consensus network
	 * @param ID the unique identifier of the resulting network
	 * @param supporting_networks the minimal number of networks that need to agree on a certain edge
	 * @param refRequired whether or not the presence of the edge in the reference network is required for it to be included in the consensus network.
	 * @param weight_cutoff the minimal value of a resulting edge for it to be included in the consensus network
	 * @param log the logger that records logging messages
	 * @param minOperator if true, the minimum of all matching edges is taken to calculate the consensus, otherwise the maximum
	 * 
	 * @return the consensus network between all input networks
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	private ConsensusNetwork calculateConsensusNetwork(Set<Network> networks, TreeEdgeOntology eo,
			NodeMapper nm, String consensus_name, int ID, int supporting_networks, boolean refRequired, double weight_cutoff, Logger log, boolean minOperator) throws IllegalArgumentException
	{
		if (networks == null || networks.isEmpty() || eo == null || nm == null || consensus_name == null)
		{
			String errormsg = "Found null parameter in calculateConsensusNetwork!";
			throw new IllegalArgumentException(errormsg);
		}
		if (supporting_networks <= 0 || supporting_networks > networks.size())
		{
			String errormsg = "The number of supporting_networks (" + supporting_networks + ") should be between 0 (excl) and the number of input networks (incl, " + networks.size() + ")";
			throw new IllegalArgumentException(errormsg);
		}
		if (mode.equals(RunMode.EDGEBYEDGE))
		{
			ConsensusNetwork on = new EdgeByEdge(log).calculateConsensusNetwork(networks, eo, nm, consensus_name, ID, supporting_networks, refRequired, weight_cutoff, minOperator);
			new NetworkCleaning(log).fullOverlapOutputCleaning(on, eo);
			return on;
		}
		System.out.println("Encountered unknown or unsupported mode: " + mode);
		return null;
	}


	/**
	 * Calculate the differential network between the reference and a set of condition-specific networks.
	 * 
	 * @param reference the reference network
	 * @param conditions a set of condition-specific networks (at least 1)
	 * @param eo the edge ontology that provides meaning to the edge types
	 * @param nm the node mapper that allows to map nodes from the one network to the other
	 * @param diff_name the name to give to the differential network
	 * @param ID the unique identifier of the resulting network
	 * @param supporting_networks the minimal number of networks that need to agree on a certain edge
	 * @param weight_cutoff the minimal value of a resulting edge for it to be included in the consensus network
	 * @param log the logger that records logging messages
	 * 
	 * @return the differential network between the two
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	private DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditions, TreeEdgeOntology eo,
			NodeMapper nm, String diff_name, int ID, int supporting_networks, double weight_cutoff, Logger log) throws IllegalArgumentException
	{
		if (reference == null || conditions == null || conditions.isEmpty() || eo == null || nm == null || diff_name == null)
		{
			String errormsg = "The edge weight should be positive!";
			throw new IllegalArgumentException(errormsg);
		}

		if (mode.equals(RunMode.EDGEBYEDGE))
		{
			DifferentialNetwork dn = new EdgeByEdge(log).calculateDiffNetwork(reference, conditions, eo, nm, diff_name, ID, supporting_networks, weight_cutoff);
			new NetworkCleaning(log).fullDifferentialOutputCleaning(dn);
			return dn;
		}
		System.out.println("Encountered unknown or unsupported mode: " + mode);
		return null;
	}
	

	/**
	 * Calculate the differential network and/or the consensus interactions, starting from a predefined runconfiguration in a project.
	 * 
	 * The calculated differential networks are added to the project directly, after cleaning previous output first.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param runID the ID of the configuration that needs to be run
	 * @param diff_name the name to give to the differential network. The consensus network will get this name + the prefix 'consensus_'
	 * @param diff_ID the unique identifier of the resulting differential network (or negative when it should not be calculated)
	 * @param consensus_ID the unique identifier of the resulting consensus network (or negative when it should not be calculated)
	 * @param weight_cutoff the minimal value of a resulting edge for it to be included in the consensus network
	 * @param minOperator if true, the minimum of all matching edges is taken to calculate the consensus, otherwise the maximum. If null, default settings will resort to min
	 * 
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateOneDifferentialNetwork(Project p, int runID, String diff_name, int diff_ID, int consensus_ID, double weight_cutoff, Boolean minOperator) throws IllegalArgumentException
	{
		TreeEdgeOntology eo = p.getEdgeOntology();
		NodeMapper nm = p.getNodeMapper();
		Logger log = p.getLogger(runID);

		RunConfiguration rc = p.getRunConfiguration(runID);
		RunOutput output = p.getOutput(runID);
		output.clean();
		
		DifferentialNetwork diff = null;
		ConsensusNetwork on = null;
		
		if (minOperator == null)
		{
			minOperator = default_MIN;
		}
		
		if (diff_ID >= 0)
		{
			if (! p.isDiffType(runID))
			{
				String errormsg = "The provided runconfiguration is not suited to calculated differential networks!";
				throw new IllegalArgumentException(errormsg);
			}
			RunDiffConfiguration drc = (RunDiffConfiguration) rc;
			ReferenceNetwork r = drc.getReferenceNetwork();
			Set<ConditionNetwork> cs = new HashSet<ConditionNetwork>(drc.getConditionNetworks());
			log.log("Calculating the differential and consensus network between " + r.getName() + " and "
					+ cs.size() + " condition-dependent network(s)");
			diff = calculateDiffNetwork(r, cs, eo, nm, diff_name, diff_ID, rc.getOverlapCutoff()-1, weight_cutoff, log);
		}
		
		if (consensus_ID >= 0)
		{
			Set<Network> inputs = new HashSet<Network>();
			inputs.addAll(rc.getInputNetworks());
			String consensus_name = consensusnameprefix + diff_name;
			log.log("Calculating the consensus network between " + inputs.size() + " input network(s)");
			on = calculateConsensusNetwork(inputs, eo, nm, consensus_name, consensus_ID, rc.getOverlapCutoff(), rc.getRefRequired(), weight_cutoff, log, minOperator);
		}
		
		if (diff != null && on != null)
		{
			output.addPair(new OutputNetworkPair(diff, on));
		}
		else if (diff != null)
		{
			output.addDifferential(diff);
		}
		else if (on != null)
		{
			output.addOverlap(on);
		}
		log.log("Done!");
	}

	/**
	 * Calculate the differential network between the reference and a set of condition-specific networks, 
	 * as well as the counterpart consensus interactions as a related ConsensusNetwork object.
	 * 
	 * The name of the differential network will be 'all_conditions_against_reference_diff'.
	 * The name of the corresponding consensus network will get an additional prefix 'consensus_' at the end.
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * The calculated differential networks and the logger object are added are added to the project directly.
	 * The logger object will first be cleaned.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param runID the ID of the configuration that needs to be run
	 * @param diff_ID the unique identifier of the resulting differential network (or negative when it should not be calculated)
	 * @param consensus_ID the unique identifier of the resulting consensus network (or negative when it should not be calculated)
	 * @param minOperator if true, the minimum of all matching edges is taken to calculate the consensus, otherwise the maximum. If null, default settings will resort to min
	 * 
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateOneDifferentialNetwork(Project p, int runID, int diff_ID, int consensus_ID, Boolean minOperator) throws IllegalArgumentException
	{
		String diff_name = diffnameprefix + "all_conditions_against_reference";
		calculateOneDifferentialNetwork(p, runID, diff_name, diff_ID, consensus_ID, default_weight_cutoff, minOperator);
	}

	/**
	 * Calculate the differential network between the reference and a set of condition-specific networks.
	 * Adds the consensus interactions as a related ConsensusNetwork object.
	 * 
	 * The name of the differential network will be 'all_conditions_against_reference_diff'.
	 * The name of the corresponding consensus network will get an additional prefix 'consensus_' at the end.
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * The calculated differential networks and the logger object are added are added to the project directly.
	 * The logger object will first be cleaned.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param runID the ID of the configuration that needs to be run
	 * @param weight_cutoff the minimal value of a resulting edge for it to be included in the consensus network
	 * @param diff_ID the unique identifier of the resulting differential network (or negative when it should not be calculated)
	 * @param consensus_ID the unique identifier of the resulting consensus network (or negative when it should not be calculated)
	 * @param minOperator if true, the minimum of all matching edges is taken to calculate the consensus, otherwise the maximum. If null, default settings will resort to min
	 * 
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateOneDifferentialNetwork(Project p, int runID, double weight_cutoff, int diff_ID, int consensus_ID, Boolean minOperator) throws IllegalArgumentException
	{
		String diff_name = diffnameprefix + "all_conditions_against_reference";
		calculateOneDifferentialNetwork(p, runID, diff_name, diff_ID, consensus_ID, weight_cutoff, minOperator);
	}

	
	/////// METHODS BETWEEN 1 REFERENCE NETWORK AND 1 CONDITION-SPECIFIC NETWORK  //////////////

	
	/**
	 * Calculate all pairwise differential networks between the reference and each condition-specific network in the project.
	 * Adds the consensus interactions as a related ConsensusNetwork object for each differential network.
	 * 
	 * The name of the differential network will be the name of the condition-specific network with prefix 'diff_'
	 * The name of the corresponding consensus network will get an additional prefix 'consensus_'.
	 * The weight cutoff will be 0.0, meaning that all edges will be included, however small their (positive) weight is.
	 * 
	 * All calculated differential networks and the logger object are added to the project directly.
	 * The logger object will first be cleaned.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param runID the ID of the configuration that needs to be run
	 * @param diffNetwork whether or not to calculate differential networks
	 * @param consensusNetwork whether or not to calculate consensus networks
	 * @param firstID the first ID that can be used for the output networks; subsequent IDs will be constructed by adding 1 each time
	 * @param minOperator if true, the minimum of all matching edges is taken to calculate the consensus, otherwise the maximum. If null, default settings will resort to min
	 * 
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateAllPairwiseDifferentialNetworks(Project p, int runID, boolean diffNetwork, boolean consensusNetwork, int firstID, Boolean minOperator) throws IllegalArgumentException
	{
		calculateAllPairwiseDifferentialNetworks(p, runID, default_weight_cutoff, diffNetwork, consensusNetwork, firstID, minOperator);
	}

	/**
	 * Calculate all pairwise differential networks and/or consensus networks between the reference and each condition-specific network in the project.
	 * The required number of supporting networks will always be 2, and will thus not be queried from the corresponding RunConfiguration object.
	 * 
	 * The name of the differential network will be the name of the condition-specific network with prefix 'diff_'
	 * The name of the consensus networks will be prefix 'consensus_' + the two names of the networks.
	 * 
	 * All calculated differential networks and the logger object are added to the project directly.
	 * The logger object and the output result will first be cleaned.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param runID the ID of the configuration that needs to be run
	 * @param weight_cutoff the minimal value of a resulting edge for it to be included in the consensus network
	 * @param diffNetwork whether or not to calculate a differential network
	 * @param consensusNetwork whether or not to calculate a consensus network
	 * @param firstID the first ID that can be used for the output networks; subsequent IDs will be constructed by adding 1 each time
	 * @param minOperator if true, the minimum of all matching edges is taken to calculate the consensus, otherwise the maximum. If null, default settings will resort to min
	 * 
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateAllPairwiseDifferentialNetworks(Project p, int runID, double weight_cutoff, boolean diffNetwork, boolean consensusNetwork, int firstID, Boolean minOperator) throws IllegalArgumentException
	{
		TreeEdgeOntology eo = p.getEdgeOntology();
		NodeMapper nm = p.getNodeMapper();
		Logger log = p.getLogger(runID);

		RunConfiguration rc = p.getRunConfiguration(runID);
		RunOutput output = p.getOutput(runID);
		output.clean();
		
		if (minOperator == null)
		{
			minOperator = default_MIN;
		}
		
		if (diffNetwork)
		{
			if (! p.isDiffType(runID))
			{
				String errormsg = "The provided runconfiguration is not suited to calculated differential networks!";
				throw new IllegalArgumentException(errormsg);
			}
			RunDiffConfiguration drc = (RunDiffConfiguration) rc;
			ReferenceNetwork r = drc.getReferenceNetwork();
			Set<ConditionNetwork> cs = new HashSet<ConditionNetwork>(drc.getConditionNetworks());
			
			for (ConditionNetwork c : cs)
			{
				String diff_name = diffnameprefix + c.getName();
				log.log("Calculating the differential and consensus network between " + r.getName() + " and " + c.getName());
				Set<ConditionNetwork> oneCs = new HashSet<ConditionNetwork>();
				oneCs.add(c);
				DifferentialNetwork diff = calculateDiffNetwork(r, oneCs, eo, nm, diff_name, firstID++, 1, weight_cutoff, log);
				
				ConsensusNetwork on = null;
				
				if (consensusNetwork)
				{
					// create a consensus name with consistent alphabetical ordering of the network names
					String consensus_name = consensusnameprefix + r.getName() + "_" + c.getName();
					
					if (c.getName().compareTo(r.getName()) < 0)
					{
						consensus_name = consensusnameprefix + c.getName() + "_" + r.getName();
					}
					
					Set<Network> inputs = new HashSet<Network>();
					inputs.add(r);
					inputs.add(c);
					on = calculateConsensusNetwork(inputs, eo, nm, consensus_name, firstID++, 2, true, weight_cutoff, log, minOperator);
				}
				
				if (on != null)
				{
					output.addPair(new OutputNetworkPair(diff, on));
				}
				else 
				{
					output.addDifferential(diff);
				}
			}
		}
		else if (consensusNetwork)
		{
			List<InputNetwork> inputs = new ArrayList<InputNetwork>(rc.getInputNetworks());
			for (int i = 0; i < inputs.size(); i++)
			{
				InputNetwork n1 = inputs.get(i);
				for (int j = i+1; j < inputs.size(); j++)
				{
					InputNetwork n2 = inputs.get(j);
					
					log.log("Calculating the consensus network between " + n1.getName() + " and " + n2.getName());
					
					// create a consensus name with consistent alphabetical ordering of the network names
					String consensus_name = consensusnameprefix + n1.getName() + "_" + n2.getName();
					if (n2.getName().compareTo(n1.getName()) < 0)
					{
						consensus_name = consensusnameprefix + n2.getName() + "_" + n1.getName();
					}
					
					Set<Network> twoInputs = new HashSet<Network>();
					twoInputs.add(n1);
					twoInputs.add(n2);
					ConsensusNetwork on = calculateConsensusNetwork(twoInputs, eo, nm, consensus_name, firstID++, 2, false, weight_cutoff, log, minOperator);
					
					output.addOverlap(on);
				}
			}
		}
		log.log("Done!");
	}
	
}
