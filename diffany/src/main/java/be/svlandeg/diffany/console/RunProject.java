package be.svlandeg.diffany.console;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.cli.CommandLine;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.io.NetworkIO;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.progress.ProgressListener;
import be.svlandeg.diffany.core.progress.StandardProgressListener;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;

/**
 * This class can run the Diffany algorithms from a {@link org.apache.commons.cli.CommandLine} object.
 * Currently, only pairwise comparisons are supported, with 1 reference network and 1 condition-dependent network.
 * 
 * @author Sofie Van Landeghem
 */
public class RunProject
{

	/**
	 * Run a Diffany analysis, depending on the parameters provided on the commandline.
	 * 
	 * @param cmd the CommandLine object containing all parameters provided by the user
	 * 
	 * @throws IllegalArgumentException when a crucial argument is missing
	 * @throws IOException when the provided directory arguments can not be read properly
	 */
	public void runAnalysis(CommandLine cmd) throws IOException, IllegalArgumentException
	{
		// TODO v.3.0: check inappropriate argument values or silently ignore and use default (as now)?

		CalculateDiff diffAlgo = new CalculateDiff();
		ProgressListener listener = null;

		Project p = new Project("Diffany-Analysis", new DefaultEdgeOntology());

		/** PARSE INPUT **/
		boolean skipHeader = readBooleanValue(cmd, DiffanyOptions.headerShort, DiffanyOptions.defaultReadHeader, "yes", "no");
		
		if (cmd.hasOption(DiffanyOptions.logShort))
		{
			listener = new StandardProgressListener(true);
		}

		File inputDir = getRequiredDir(cmd, DiffanyOptions.inputShort);
		ReferenceNetwork refNet = null;
		Set<ConditionNetwork> conditionNets = new HashSet<ConditionNetwork>();
		 
		Set<InputNetwork> inputnetworks = NetworkIO.readGenericInputNetworksFromSubdirs(inputDir, skipHeader, false);
		
		if (inputnetworks == null || inputnetworks.size() < 2)
		{
			String msg = "Could not read all required input networks from " + inputDir;
			throw new IllegalArgumentException(msg);
		}
		
		for (InputNetwork net : inputnetworks)
		{
			if (net instanceof ReferenceNetwork)
			{
				if (refNet != null)
				{
					String msg = "Found more than 1 reference network at " + inputDir;
					throw new IllegalArgumentException(msg);
				}
				refNet = (ReferenceNetwork) net;
			}
			else if (net instanceof ConditionNetwork)
			{
				conditionNets.add((ConditionNetwork) net);
			}
			else
			{
				String msg = "Found a strange input network: " + net;
				throw new IllegalArgumentException(msg);
			}
		}

		boolean runDiff = readBooleanValue(cmd, DiffanyOptions.runDiff, DiffanyOptions.defaultRunDiff, "yes", "no");
		boolean runCons = readBooleanValue(cmd, DiffanyOptions.runCons, DiffanyOptions.defaultRunCons, "yes", "no");
		
		if (runDiff && refNet == null)
		{
			String msg = "Could not read the reference network at " + inputDir;
			throw new IllegalArgumentException(msg);
		}

		int nextID = inferNextID(cmd, inputnetworks);

		Double cutoff = null;
		if (cmd.hasOption(DiffanyOptions.cutoffShort))
		{
			cutoff = Double.parseDouble(cmd.getOptionValue(DiffanyOptions.cutoffShort));
		}

		boolean minOperator = readBooleanValue(cmd, DiffanyOptions.operatorShort, DiffanyOptions.defaultMinOperator, "min", "max");
		boolean modePairwise = readBooleanValue(cmd, DiffanyOptions.modeShort, DiffanyOptions.defaultModePairwise, "pairwise", "all");

		/** THE ACTUAL ALGORITHM **/
		boolean cleanInput = true;
		Integer runID = p.addRunConfiguration(refNet, conditionNets, cleanInput, listener);

		if (modePairwise)
		{
			diffAlgo.calculateAllPairwiseDifferentialNetworks(p, runID, cutoff, runDiff, runCons, nextID, minOperator, listener);
		}
		else
		{
			int diffID = -1;
			if (runDiff)
			{
				diffID = nextID++;
			}
			int consensusID = -1;
			if (runCons)
			{
				consensusID = nextID++;
			}

			// default names for the output networks will be generated
			diffAlgo.calculateOneDifferentialNetwork(p, runID, cutoff, null, null, diffID, consensusID, minOperator, listener);
		}

		/** WRITE NETWORK OUTPUT **/
		RunOutput output = p.getOutput(runID);
		boolean writeHeaders = true;
		File outputDir = getRequiredDir(cmd, DiffanyOptions.outputShort);

		for (DifferentialNetwork diffNet : output.getDifferentialNetworks())
		{
			File diffDir = new File(outputDir, "Differential_network_" + diffNet.getID());
			NetworkIO.writeNetworkToDir(diffNet, diffDir, writeHeaders);
		}

		for (ConsensusNetwork consensusNet : output.getConsensusNetworks())
		{
			File consensusDir = new File(outputDir, "Consensus_network_" + consensusNet.getID());
			NetworkIO.writeNetworkToDir(consensusNet, consensusDir, writeHeaders);
		}
	}

	/**
	 * Parse the required first ID for the output network from the optional arguments,
	 * or determine it from the given input networks.
	 */
	private int inferNextID(CommandLine cmd, Set<InputNetwork> inputNetworks)
	{
		int nextID = -1;
		if (cmd.hasOption(DiffanyOptions.nextID))
		{
			nextID = Integer.parseInt(cmd.getOptionValue(DiffanyOptions.nextID));
		}
		else
		{
			for (InputNetwork input : inputNetworks)
			{
				nextID = Math.max(nextID, input.getID() + 1);
			}
		}
		return nextID;
	}

	/**
	 * Retrieve a File object representing a directory given by a value on the command line
	 * 
	 * @param cmd the parsed command line arguments
	 * @param key the key used to specify the directory on the command line
	 * 
	 * @return the value of the key, represented as a directory object
	 * @throws IllegalArgumentException when a crucial argument is missing
	 */
	private File getRequiredDir(CommandLine cmd, String key) throws IllegalArgumentException
	{
		String inputDir = cmd.getOptionValue(key);
		if (inputDir == null)
		{
			throw new IllegalArgumentException("Fatal error: please provide a valid directory pointer for " + key);
		}
		return new File(inputDir);
	}
	
	/**
	 * Retrieve a boolean value for a given key 
	 * 
	 * @param cmd the parsed command line arguments
	 * @param key the key used to specify the boolean value
	 * 
	 * @return the value of the key, represented as a boolean
	 */
	private Boolean readBooleanValue(CommandLine cmd, String key, boolean defaultValue, String trueValue, String falseValue)
	{
		boolean result = defaultValue;
		if (cmd.hasOption(key))
		{
			String value = cmd.getOptionValue(key);
			if (value != null && value.trim().equals(trueValue))
			{
				result = true;
			}
			if (value != null && value.trim().equals(falseValue))
			{
				result = false;
			}
		}
		return result;
	}

}
