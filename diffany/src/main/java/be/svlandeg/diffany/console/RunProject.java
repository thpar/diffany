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
		CalculateDiff diffAlgo = new CalculateDiff();
		ProgressListener listener = null;

		Project p = new Project("Diffany-Analysis", new DefaultEdgeOntology());

		/** PARSE INPUT **/
		boolean skipHeader = true;
		if (cmd.hasOption(DiffanyOptions.logShort))
		{
			listener = new StandardProgressListener(true);
		}
		
		File refDir = getRequiredDir(cmd, DiffanyOptions.refShort);
		ReferenceNetwork refNet = NetworkIO.readReferenceNetworkFromDir(refDir, skipHeader);

		File condDir = getRequiredDir(cmd, DiffanyOptions.condShort);
		ConditionNetwork condNet = NetworkIO.readConditionNetworkFromDir(condDir, skipHeader);
		
		Set<InputNetwork> inputnetworks = new HashSet<InputNetwork>();
		inputnetworks.add(refNet);
		inputnetworks.add(condNet);

		/* it's no problem if these are null, CalculateDiff will then resort to a default option */
		String diffname = cmd.getOptionValue(DiffanyOptions.diffNameShort);	
		String consensusname = cmd.getOptionValue(DiffanyOptions.consNameShort);	
		
		int diffID = inferDiffID(cmd, inputnetworks);
		int consensusID = inferConsID(cmd, diffID);

		Double cutoff = null;
		if (cmd.hasOption(DiffanyOptions.cutoffShort))
		{
			cutoff = Double.parseDouble(cmd.getOptionValue(DiffanyOptions.cutoffShort));
		}
		
		boolean minOperator = DiffanyOptions.defaultMinOperator;
		if (cmd.hasOption(DiffanyOptions.operatorShort))
		{
			String operator = cmd.getOptionValue(DiffanyOptions.operatorShort);
			if (operator != null && operator.trim().equals("max"))
			{
				minOperator = false;
			}
			if (operator != null && operator.trim().equals("min"))
			{
				minOperator = true;
			}
			// TODO v.3.0: check inappropriate argument value or silently ignore and use default (as now)?
		}
		
		boolean modePairwise = DiffanyOptions.defaultModePairwise;
		if (cmd.hasOption(DiffanyOptions.modeShort))
		{
			String mode = cmd.getOptionValue(DiffanyOptions.modeShort);
			if (mode != null && mode.trim().equals("pairwise"))
			{
				modePairwise = true;
			}
			if (mode != null && mode.trim().equals("all"))
			{
				modePairwise = false;
			}
			// TODO v.3.0: check inappropriate argument value or silently ignore and use default (as now)?
		}

		/** THE ACTUAL ALGORITHM **/
		boolean cleanInput = true;
		Integer runID = p.addRunConfiguration(refNet, condNet, cleanInput, listener);

		// TODO v2.1: allow to change mode pairwise vs. differential
		if (modePairwise)
		{
			// TODO
		}
		else
		{
			diffAlgo.calculateOneDifferentialNetwork(p, runID, cutoff, diffname, consensusname, diffID, consensusID, minOperator, listener);
		}

		/** WRITE NETWORK OUTPUT **/
		RunOutput output = p.getOutput(runID);
		boolean writeHeaders = true;
		File outputDir = getRequiredDir(cmd, DiffanyOptions.outputShort);
		
		for (DifferentialNetwork diffNet : output.getDifferentialNetworks())
		{			
			File diffDir = new File (outputDir, "Reference_network_" + diffNet.getID());
			NetworkIO.writeNetworkToDir(diffNet, diffDir, writeHeaders);
		}
		
		for (ConsensusNetwork consensusNet : output.getConsensusNetworks())
		{
			File consensusDir = new File (outputDir, "Consensus_network_" + consensusNet.getID());
			NetworkIO.writeNetworkToDir(consensusNet, consensusDir, writeHeaders);
		}
	}

	/**
	 * Parse the required ID for the differential network from the optional arguments,
	 * or determine it from the given input networks.
	 */
	private int inferDiffID(CommandLine cmd, Set<InputNetwork> inputNetworks)
	{
		int diffID = -1;
		if (cmd.hasOption(DiffanyOptions.diffID))
		{
			diffID = Integer.parseInt(cmd.getOptionValue(DiffanyOptions.diffID));
		}
		else
		{
			for (InputNetwork input : inputNetworks)
			{
				diffID = Math.max(diffID, input.getID()+1);
			}
		}
		return diffID;
	}
	
	/**
	 * Parse the required ID for the consensus network from the optional arguments,
	 * or determine it from the given ID for the differential network.
	 */
	private int inferConsID(CommandLine cmd, Integer diffID)
	{
		int consID = -1;
		if (cmd.hasOption(DiffanyOptions.consID))
		{
			consID = Integer.parseInt(cmd.getOptionValue(DiffanyOptions.consID));
		}
		else
		{
			consID = diffID+1;
		}
		return consID;
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
		String refDir = cmd.getOptionValue(key);
		if (refDir == null)
		{
			throw new IllegalArgumentException("Fatal error: please provide a valid " + DiffanyOptions.refShort + " directory pointer!");
		}
		return new File(refDir);
	}

}
