package be.svlandeg.diffany.core.algorithms;

/*
 * #%L
 * Diffany
 * %%
 * Copyright (C) 2014 PSB/UGent - Sofie Van Landeghem and Thomas Van Parys
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */

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
import be.svlandeg.diffany.core.progress.ProgressListener;
import be.svlandeg.diffany.core.progress.ScheduledTask;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunConfiguration;
import be.svlandeg.diffany.core.project.RunDiffConfiguration;
import be.svlandeg.diffany.core.project.RunOutput;
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

	protected double default_weight_cutoff = 0.0;
	
	protected static final String diffname_default_all = "diff_all_conditions_against_reference";
	protected static String consensusname_default_all = "consensus_all_input_networks";
	
	protected static final String diffname_pairwise_prefix = "diff_";
	protected static final String consensusname_pairwise_prefix = "consensus_";

	/**
	 * A run mode determines how differential networks will be calculated.
	 * Currently, only the edge by edge option is supported, but this may be extended in the future.
	 */
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
	 * @param consensus_name the name to give to the consensus network (can be null)
	 * @param ID the unique identifier of the resulting network
	 * @param supportingCutoff the minimal number of networks that need to agree on a certain edge
	 * @param refRequired whether or not the presence of the edge in the reference network is required for it to be included in the consensus network.
	 * @param weightCutoff the minimal value of a resulting edge for it to be included in the differential or consensus network 
	 * @param log the logger that records logging messages
	 * @param minOperator if true, the minimum of all matching edges is taken to calculate the consensus, otherwise the maximum
	 * @param task the task object that keeps track of the progress of this calculation (can be null)
	 * 
	 * @return the consensus network between all input networks
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	private ConsensusNetwork calculateConsensusNetwork(Set<Network> networks, TreeEdgeOntology eo, String consensus_name, int ID, 
			int supportingCutoff, boolean refRequired, double weightCutoff, Logger log, boolean minOperator, ScheduledTask task) throws IllegalArgumentException
	{
		if (networks == null || networks.isEmpty() || eo == null)
		{
			String errormsg = "Found null parameter in calculateConsensusNetwork!";
			errormsg += " (null values: networks=" + (networks == null) + " / eo=" + (eo == null) 
					+ " / consensus_name="  + (consensus_name == null) + " / networks empty="  + (networks.isEmpty());
			throw new IllegalArgumentException(errormsg);
		}
		if (consensus_name == null)
		{
			consensus_name = consensusname_default_all;
		}
		if (supportingCutoff <= 0 || supportingCutoff > networks.size())
		{
			String errormsg = "The number of supportingCutoff (" + supportingCutoff + ") should be between 0 (excl) and the number of input networks (incl, " + networks.size() + ")";
			throw new IllegalArgumentException(errormsg);
		}
		ScheduledTask cleaningTask = null;
		ScheduledTask calculatingTask = null;
		if (task != null)
		{
			int totalTicks = task.ticksToGo();
			
			int cleanTicks = totalTicks / 10;
			int calculationTicks = totalTicks - cleanTicks;
			
			cleaningTask = new ScheduledTask(task.getListener(), cleanTicks);
			calculatingTask = new ScheduledTask(task.getListener(), calculationTicks);
		}
		if (mode.equals(RunMode.EDGEBYEDGE))
		{
			ConsensusNetwork on = new EdgeByEdge(log).calculateConsensusNetwork(networks, eo, consensus_name, ID, supportingCutoff, refRequired, weightCutoff, minOperator, calculatingTask);
			new NetworkCleaning(log).fullConsensusOutputCleaning(on, eo, cleaningTask);
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
	 * @param diff_name the name to give to the differential network (can be null)
	 * @param ID the unique identifier of the resulting network
	 * @param supportingCutoff the minimal number of networks that need to agree on a certain edge
	 * @param weightCutoff the minimal value of a resulting edge for it to be included in the differential or consensus network
	 * @param log the logger that records logging messages
	 * @param task the task object that keeps track of the progress of this calculation (can be null)
	 * 
	 * @return the differential network between the two
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	private DifferentialNetwork calculateDiffNetwork(ReferenceNetwork reference, Set<ConditionNetwork> conditions, TreeEdgeOntology eo,
			String diff_name, int ID, int supportingCutoff, double weightCutoff, Logger log, ScheduledTask task) throws IllegalArgumentException
	{
		if (reference == null || conditions == null || conditions.isEmpty() || eo == null)
		{
			String errormsg = "Found null parameter in calculateDiffNetwork!";
			errormsg += " (null values: reference=" + (reference == null) + " / conditions="  + (conditions == null) 
					+ " / eo=" + (eo == null) + " / conditions empty="  + (conditions.isEmpty());
			throw new IllegalArgumentException(errormsg);
		}
		if (diff_name == null)
		{
			diff_name = diffname_default_all;
		}
		ScheduledTask cleaningTask = null;
		ScheduledTask calculatingTask = null;
		if (task != null)
		{
			int totalTicks = task.ticksToGo();
			int cleanTicks = totalTicks / 10;
			int calculationTicks = totalTicks - cleanTicks;
			
			cleaningTask = new ScheduledTask(task.getListener(), cleanTicks);
			calculatingTask = new ScheduledTask(task.getListener(), calculationTicks);
		}
		if (mode.equals(RunMode.EDGEBYEDGE))
		{
			DifferentialNetwork dn = new EdgeByEdge(log).calculateDiffNetwork(reference, conditions, eo, diff_name, ID, supportingCutoff, weightCutoff, calculatingTask);
			new NetworkCleaning(log).fullDifferentialOutputCleaning(dn, eo, cleaningTask);
			return dn;
		}
		System.out.println("Encountered unknown or unsupported mode: " + mode);
		return null;
	}
	

	/**
	 * Calculate the differential network and/or the consensus network, starting from a predefined runconfiguration in a project.
	 * 
	 * The calculated differential networks are added to the project directly, after cleaning previous output first.
	 * 
	 * @param p the project which stores the reference and condition-specific networks
	 * @param runID the ID of the configuration that needs to be run
	 * @param diff_name the name to give to the differential network (can be null, then a default name will be constructed)
	 * @param consensus_name the name to give to the differential network (can be null, then a default name will be constructed)
	 * @param diff_ID the unique identifier of the resulting differential network (or negative when it should not be calculated)
	 * @param consensus_ID the unique identifier of the resulting consensus network (or negative when it should not be calculated)
	 * @param weight_cutoff the minimal value of a resulting edge for it to be included in the differential or consensus network (can be null, then a default value will be used)
	 * @param minOperator if true, the minimum of all matching edges is taken to calculate the consensus, otherwise the maximum (can be null, then 'min' will be used by default)
	 * @param progressListener the listener that will be updated about the progress of this calculation (can be null)
	 * 
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateOneDifferentialNetwork(Project p, int runID, Double weight_cutoff, String diff_name, String consensus_name, int diff_ID, int consensus_ID, 
			Boolean minOperator, ProgressListener progressListener) throws IllegalArgumentException
	{
		if (weight_cutoff == null)
		{
			weight_cutoff = default_weight_cutoff;
		}
		
		TreeEdgeOntology eo = p.getEdgeOntology();
		Logger log = p.getLogger(runID);

		RunConfiguration rc = p.getRunConfiguration(runID);
		RunOutput output = p.getOutput(runID);
		output.clean();
		
		DifferentialNetwork diff = null;
		ConsensusNetwork cn = null;
		
		if (minOperator == null)
		{
			minOperator = default_MIN;
		}
		
		int totalTicks = 100;
		ScheduledTask diffTask = null;
		ScheduledTask consensusTask = null;
		if (progressListener != null)
		{
			progressListener.reset(totalTicks);
			
			int diffTicks = totalTicks;
			int consensusTicks = totalTicks;
			
			if (diff_ID >= 0)
			{
				consensusTicks = consensusTicks / 2;
			}
			if (consensus_ID >= 0)
			{
				diffTicks = totalTicks - consensusTicks;
			}
			diffTask = new ScheduledTask(progressListener, diffTicks);
			consensusTask = new ScheduledTask(progressListener, consensusTicks);
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
			log.log("Calculating the 1-all differential network between " + r.getName() + " and "
					+ cs.size() + " condition-dependent network(s) - settings: support cutoff = " + rc.getSupportCutoff() + ", weight cutoff = " + weight_cutoff);
			diff = calculateDiffNetwork(r, cs, eo, diff_name, diff_ID, rc.getSupportCutoff()-1, weight_cutoff, log, diffTask);
		}
		
		if (consensus_ID >= 0)
		{
			Set<Network> inputs = new HashSet<Network>();
			inputs.addAll(rc.getInputNetworks());
			log.log("Calculating the 1-all consensus network between " + inputs.size() + " input network(s) - settings: support cutoff = " + rc.getSupportCutoff() + ", weight cutoff = " + weight_cutoff + ", reference required = " + rc.getRefRequired() + ", minOperator = " + minOperator);
			cn = calculateConsensusNetwork(inputs, eo, consensus_name, consensus_ID, rc.getSupportCutoff(), rc.getRefRequired(), weight_cutoff, log, minOperator, consensusTask);
		}
		
		if (diff != null && cn != null)
		{
			output.addPair(new OutputNetworkPair(diff, cn));
		}
		else if (diff != null)
		{
			output.addDifferential(diff);
		}
		else if (cn != null)
		{
			output.addConsensus(cn);
		}
		log.log("Done!");
	}

	
	/////// METHODS BETWEEN 1 REFERENCE NETWORK AND 1 CONDITION-SPECIFIC NETWORK  //////////////


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
	 * @param weightCutoff the minimal value of a resulting edge for it to be included in the differential or consensus network (can be null)
	 * @param diffNetwork whether or not to calculate a differential network
	 * @param consensusNetwork whether or not to calculate a consensus network
	 * @param firstID the first ID that can be used for the output networks; subsequent IDs will be constructed by adding 1 each time
	 * @param minOperator if true, the minimum of all matching edges is taken to calculate the consensus, otherwise the maximum. If null, default settings will resort to min
	 * @param progressListener the listener that will be updated about the progress of this calculation (can be null)
	 * 
	 * @throws IllegalArgumentException if any of the crucial fields in the project are null
	 */
	public void calculateAllPairwiseDifferentialNetworks(Project p, int runID, Double weightCutoff, boolean diffNetwork, boolean consensusNetwork, int firstID, Boolean minOperator, ProgressListener progressListener) throws IllegalArgumentException
	{
		if (weightCutoff == null)
		{
			weightCutoff = default_weight_cutoff;
		}
		
		TreeEdgeOntology eo = p.getEdgeOntology();
		Logger log = p.getLogger(runID);

		RunConfiguration rc = p.getRunConfiguration(runID);
		RunOutput output = p.getOutput(runID);
		output.clean();
		
		int tasks = 0;
		int conditions = rc.getInputNetworks().size();
		if (diffNetwork)
		{
			conditions = conditions - 1;
			tasks += conditions ;			// as many pairwise comparisons as there are conditions
			if (consensusNetwork)
			{
				tasks += conditions ;		// and equally many consensus networks
			}
		}
		else if (consensusNetwork)
		{
			tasks = (conditions * (conditions-1)) / 2;		// all possible pairwise combinations of the conditions
		}
		
		int ticksPerTask = 100;
		if (progressListener != null)
		{
			progressListener.reset(tasks * ticksPerTask);
		}
		
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
				ScheduledTask diffTask = null;
				if (progressListener != null)
				{
					diffTask = new ScheduledTask(progressListener, ticksPerTask);
				}
				
				String diff_name = diffname_pairwise_prefix + c.getName();
				log.log("Calculating the pairwise differential network between " + r.getName() + " and " + c.getName() + " - settings: weight cutoff = " + weightCutoff);
				Set<ConditionNetwork> oneCs = new HashSet<ConditionNetwork>();
				oneCs.add(c);
				DifferentialNetwork diff = calculateDiffNetwork(r, oneCs, eo, diff_name, firstID++, 1, weightCutoff, log, diffTask);
				
				ConsensusNetwork on = null;
				
				if (consensusNetwork)
				{
					// create a consensus name with consistent alphabetical ordering of the network names
					String consensus_name = consensusname_pairwise_prefix + r.getName() + "_" + c.getName();
					
					if (c.getName().compareTo(r.getName()) < 0)
					{
						consensus_name = consensusname_pairwise_prefix + c.getName() + "_" + r.getName();
					}
					
					Set<Network> inputs = new HashSet<Network>();
					inputs.add(r);
					inputs.add(c);
					
					ScheduledTask consTask = null;
					if (progressListener != null)
					{
						consTask = new ScheduledTask(progressListener, ticksPerTask);
					}
					log.log("Calculating the pairwise consensus network between " + r.getName() + " and " + c.getName() + " - settings: weight cutoff = " + weightCutoff + ", minOperator = " + minOperator);
					
					on = calculateConsensusNetwork(inputs, eo, consensus_name, firstID++, 2, true, weightCutoff, log, minOperator, consTask);
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
					ScheduledTask consTask = null;
					if (progressListener != null)
					{
						consTask = new ScheduledTask(progressListener, ticksPerTask);
					}
					
					InputNetwork n2 = inputs.get(j);
					
					log.log("Calculating the consensus network between " + n1.getName() + " and " + n2.getName() + " - settings: weight cutoff = " + weightCutoff + ", minOperator = " + minOperator);
					
					// create a consensus name with consistent alphabetical ordering of the network names
					String consensus_name = consensusname_pairwise_prefix + n1.getName() + "_" + n2.getName();
					if (n2.getName().compareTo(n1.getName()) < 0)
					{
						consensus_name = consensusname_pairwise_prefix + n2.getName() + "_" + n1.getName();
					}
					
					Set<Network> twoInputs = new HashSet<Network>();
					twoInputs.add(n1);
					twoInputs.add(n2);
					ConsensusNetwork cn = calculateConsensusNetwork(twoInputs, eo, consensus_name, firstID++, 2, false, weightCutoff, log, minOperator, consTask);
					
					output.addConsensus(cn);
				}
			}
		}
		log.log("Done!");
	}
	
}
