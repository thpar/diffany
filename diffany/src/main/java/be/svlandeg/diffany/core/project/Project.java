package be.svlandeg.diffany.core.project;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import be.svlandeg.diffany.core.algorithms.NetworkCleaning;
import be.svlandeg.diffany.core.algorithms.Unification;
import be.svlandeg.diffany.core.io.ProjectIO;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.progress.ProgressListener;
import be.svlandeg.diffany.core.progress.ScheduledTask;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

/**
 * A project consists of a number of Diffany runs and the ontology settings within one user session.
 * 
 * Additionally, a project links to an {@link EdgeOntology} that defines the semantics of edge types.
 * 
 * Project data can be saved and loaded through the {@link ProjectIO} class.
 * 
 * @author Sofie Van Landeghem
 */
public class Project
{
	
	protected String name;

	protected TreeEdgeOntology edgeOntology;
	
	// List of runs; their IDs are given by the index in this list
	protected List<Run> runs;
	
	
	/**
	 * Create a new project with a default node mapper and a default edge ontology that can interpret the differential edges.
	 * 
	 * @param name the name of this project (not null!)
	 * 
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public Project(String name) throws IllegalArgumentException
	{
		this(name, new DefaultEdgeOntology());
	}

	/**
	 * Create a new project with a node mapper and an edge ontology that can interpret the differential edges.
	 * 
	 * @param name the name of this project (not null!)
	 * @param edgeOntology the edge ontology (not null!)
	 * 
	 * @throws IllegalArgumentException if any of the restrictions above are not fulfilled
	 */
	public Project(String name, TreeEdgeOntology edgeOntology) throws IllegalArgumentException
	{
		if (name == null)
		{
			String errormsg = "The name of a project should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.name = name;
		
		setEdgeOntology(edgeOntology);
		
		runs = new ArrayList<Run>();
	}
	
	/**
	 * Retrieve the number of runs in this project
	 * @return the number of run IDs
	 */
	public Integer getNumberOfRuns()
	{
		return runs.size();
	}
	
	/**
	 * Get a specific RunConfiguration by its unique run ID
	 * 
	 * @param runID the (unique) run ID of the configuration
	 * @return the corresponding RunConfiguration
	 * @throws IllegalArgumentException if the run ID is invalid (i.e. non-existing)
	 */
	public RunConfiguration getRunConfiguration(int runID) throws IllegalArgumentException
	{
		if (runID < 0 || runID >= runs.size())
		{
			String errormsg = "Unknown run ID " + runID + " in Project " + name;
			throw new IllegalArgumentException(errormsg);
		}
		return runs.get(runID).configuration;
	}
	
	
	/**
	 * Add a new RunConfiguration to this project. 
	 * There is NO check whether or not this configuration was added previously, it will simply be duplicated with a new ID!
	 * 
	 * @param reference the reference network (not null!)
	 * @param condition the condition-specific networks (not null!)
	 * @param cleanInput whether or not to clean the input networks
	 * @param progressListener the listener that will be updated about the progress of this calculation (can be null)
	 * 
	 * @return the unique run ID assigned to the new RunConfiguration in this project
	 */
	public int addRunConfiguration(ReferenceNetwork reference, ConditionNetwork condition, boolean cleanInput, ProgressListener progressListener)
	{
		Set<ConditionNetwork> cs = new HashSet<ConditionNetwork>();
		cs.add(condition);
		return addRunConfiguration(reference, cs, 2, cleanInput, progressListener);
	}
	
	/**
	 * Add a new RunConfiguration to this project, automatically registering the networks to this project.
	 * This configuration will allow the calculation of both differential and consensus networks.
	 * There is NO check whether or not this configuration was added previously, it will simply be duplicated with a new ID!
	 * This option will use the default requirement of all networks having a consensus edge before it can be include in the ConsensusNetwork.
	 * 
	 * @param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * @param cleanInput whether or not to clean the input networks
	 * @param progressListener the listener that will be updated about the progress of this calculation (can be null)
	 * 
	 * @return the unique run ID assigned to the new RunConfiguration in this project
	 */
	public int addRunConfiguration(ReferenceNetwork reference, Set<ConditionNetwork> conditions, boolean cleanInput, ProgressListener progressListener)
	{
		return addRunConfiguration(reference, conditions, conditions.size() + 1, cleanInput, progressListener);
	}
	
	/**
	 * Add a new RunConfiguration to this project, automatically registering the networks to this project.
	 * This configuration will allow the calculation of both differential and consensus networks.
	 * There is NO check whether or not this configuration was added previously, it will simply be duplicated with a new ID!
	 * 
	 * @param reference the reference network (not null!)
	 * @param conditions the condition-specific networks (at least 1!)
	 * @param supportingCutoff the number of networks that should at least match for consensus to be defined: min. 2, max conditions.size + 1.
	 * @param cleanInput whether or not to clean the input networks
	 * @param progressListener the listener that will be updated about the progress of this calculation (can be null)
	 * 
	 * @return the unique run ID assigned to the new RunConfiguration in this project
	 * @throws IllegalArgumentException when the IDs of the provided networks are not unique
	 */
	public int addRunConfiguration(ReferenceNetwork reference, Set<ConditionNetwork> conditions, int supportingCutoff, boolean cleanInput, ProgressListener progressListener)
	{
		Logger logger = new Logger();
		logger.log("Analysing the reference and condition-specific network(s) ");
		
		RunConfiguration rc = null;
		
		if (cleanInput)
		{
			/* It is necessary to first register before cleaning, otherwise some interaction types may be unknown by the NetworkCleaning object */
			int tasks = conditions.size() + 1;
			int ticksPerTask = 100;
			ScheduledTask diffTask = null;
			if (progressListener != null)
			{
				progressListener.reset(tasks * ticksPerTask);
				diffTask = new ScheduledTask(progressListener, ticksPerTask);
			}
			registerSourceNetwork(reference, logger);
			ReferenceNetwork cleanRef = new NetworkCleaning(logger).fullInputRefCleaning(reference, edgeOntology, diffTask);
			
			Set<ConditionNetwork> cleanConditions = new HashSet<ConditionNetwork>();
			for (ConditionNetwork conNet : conditions)
			{
				ScheduledTask consTask = null;
				if (progressListener != null)
				{
					consTask = new ScheduledTask(progressListener, ticksPerTask);
				}
				registerSourceNetwork(conNet, logger);
				ConditionNetwork cleanCon = new NetworkCleaning(logger).fullInputConditionCleaning(conNet, edgeOntology, consTask);
				cleanConditions.add(cleanCon);
			}
			rc = new RunDiffConfiguration(cleanRef, cleanConditions, supportingCutoff);
		}
		else
		{
			registerSourceNetwork(reference, logger);
			for (ConditionNetwork conNet : conditions)
			{
				registerSourceNetwork(conNet, logger);
			}
			rc = new RunDiffConfiguration(reference, conditions, supportingCutoff);
		}
		
		int nextID = runs.size();
		Run run = new Run(this, nextID, rc, true, logger);
		
		boolean IDsOK = run.checkInputIDs();
		if (! IDsOK)
		{
			String errormsg = "The IDs of the input networks should be unique!";
			throw new IllegalArgumentException(errormsg);
		}
		
		runs.add(run);
		
		return nextID;
	}
	
	/**
	 * Add a new RunConfiguration to this project, automatically registering the networks to this project.
	 * This configuration will NOT be able to calculate differential networks; only consensus networks.
	 * There is NO check whether or not this configuration was added previously, it will simply be duplicated with a new ID!
	 * This option will use the default requirement of all networks having a consensus edge before it can be include in the ConsensusNetwork.
	 * 
	 * @param inputNetworks all the input networks (at least 1!)
	 * @param progressListener the listener that will be updated about the progress of this calculation (can be null)
	 * 
	 * @return the unique run ID assigned to the new RunConfiguration in this project
	 */
	public int addRunConfiguration(Set<InputNetwork> inputNetworks, ProgressListener progressListener)
	{
		return addRunConfiguration(inputNetworks, inputNetworks.size(), false, progressListener);
	}
	
	/**
	 * Add a new RunConfiguration to this project, automatically registering the networks to this project.
	 * This configuration will NOT be able to calculate differential networks; only consensus networks.
	 * There is NO check whether or not this configuration was added previously, it will simply be duplicated with a new ID!
	 * 
	 * @param inputNetworks all the input networks (at least 1!)
	 * @param supportingCutoff the required number of input networks that need to match for a consensus edge to be present
	 * @param refRequired whether or not the presence of the edge in the reference network is required for it to be included in the consensus network
	 * @param progressListener the listener that will be updated about the progress of this calculation (can be null)
	 * 
	 * @return the unique run ID assigned to the new RunConfiguration in this project
	 */
	public int addRunConfiguration(Set<InputNetwork> inputNetworks, int supportingCutoff, boolean refRequired, ProgressListener progressListener)
	{
		Logger logger = new Logger();
		logger.log("Analysing the input networks ");
		
		int tasks = inputNetworks.size();
		int ticksPerTask = 100;
		if (progressListener != null)
		{
			progressListener.reset(tasks * ticksPerTask);
		}
		
		Set<InputNetwork> cleanNetworks = new HashSet<InputNetwork>();
		for (InputNetwork inputNet : inputNetworks)
		{
			ScheduledTask task = null;
			if (progressListener != null)
			{
				task = new ScheduledTask(progressListener, ticksPerTask);
			}
			
			registerSourceNetwork(inputNet, logger);
			InputNetwork cleanNet = new NetworkCleaning(logger).fullInputCleaning(inputNet, edgeOntology, task);
			cleanNetworks.add(cleanNet);
		}
		
		RunConfiguration rc = new RunConfiguration(cleanNetworks, supportingCutoff, refRequired);
		
		int nextID = runs.size();
		Run run = new Run(this, nextID, rc, false, logger);
		
		boolean IDsOK = run.checkInputIDs();
		if (! IDsOK)
		{
			String errormsg = "The IDs of the input networks should be unique!";
			throw new IllegalArgumentException(errormsg);
		}
		
		runs.add(run);
		
		return nextID;
	}
	
	/**
	 * Get the output for a specific run ID. The RunOutput object will be empty if the corresponding run configuration was not deployed.
	 * @param runID the run ID
	 * @return the output of a specific run
	 */
	public RunOutput getOutput(int runID)
	{
		return runs.get(runID).output;
	}
	
	/**
	 * Get the logger for a specific run ID. The logs will be empty if the corresponding run configuration was not deployed.
	 * @param runID the run ID
	 * @return the logs of a specific run.
	 */
	public Logger getLogger(int runID)
	{
		return runs.get(runID).logger;
	}
	
	/**
	 * Get the type of the run configuration by ID: true if it can calculate differential networks, false otherwise (only consensus).
	 * @param runID the run ID
	 * @return the type of a specific run configuration.
	 */
	public boolean isDiffType(int runID)
	{
		return runs.get(runID).type;
	}
	
	
	/**
	 * Register a new source network to this project, updating the EdgeOntology and NodeMapper instances.
	 * A source network is a reference, condition-specific, or a consensus network.
	 * There is NO check whether or not this network was added previously.
	 * 
	 * @param source the new Network that will be used within this project
	 * @param logger the logger instance that can store log messages
	 */
	public void registerSourceNetwork(Network source, Logger logger)
	{
		new Unification(logger).expandEdgeOntology(source.getEdges(), edgeOntology);
	}

	/**
	 * Return the name of this project
	 * @return the name of this project
	 */
	public String getName()
	{
		return name;
	}


	/**
	 * Set the edge ontology for this project.
	 * (currently not a public method - changes to it would influence the differential networks (TODO v3.0))
	 * 
	 * @param edgeOntology the edge ontology (not null!)
	 * @throws IllegalArgumentException if the edgeOntology is null
	 */
	private void setEdgeOntology(TreeEdgeOntology edgeOntology) throws IllegalArgumentException
	{
		if (edgeOntology == null)
		{
			String errormsg = "The edge ontology should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.edgeOntology = edgeOntology;

	}
	
	/**
	 * Get the edge ontology of this project, which can translate edge types to categories and assign semantics to the categories.
	 * 
	 * @return the edge ontology used in this project (should not be null)
	 */
	public TreeEdgeOntology getEdgeOntology()
	{
		return edgeOntology;
	}


}
