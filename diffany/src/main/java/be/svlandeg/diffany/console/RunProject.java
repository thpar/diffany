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
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.progress.ProgressListener;
import be.svlandeg.diffany.core.progress.StandardProgressListener;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.core.project.LogEntry;
import be.svlandeg.diffany.core.project.Logger;
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
		ProgressListener listener = new StandardProgressListener(true);

		Project p = new Project("Diffany-Analysis", new DefaultEdgeOntology());

		/** PARSE INPUT **/
		boolean toLog = false;
		boolean skipHeader = true;
		if (cmd.hasOption(DiffanyOptions.logShort))
		{
			toLog = true;
		}
		
		File refDir = getRequiredDir(cmd, DiffanyOptions.refShort);
		ReferenceNetwork refNet = NetworkIO.readReferenceNetworkFromDir(refDir, skipHeader);

		File condDir = getRequiredDir(cmd, DiffanyOptions.conShort);
		ConditionNetwork condNet = NetworkIO.readConditionNetworkFromDir(condDir, skipHeader);
		
		Set<InputNetwork> inputnetworks = new HashSet<InputNetwork>();
		inputnetworks.add(refNet);
		inputnetworks.add(condNet);

		String name = cmd.getOptionValue(DiffanyOptions.diffnameShort);
		int diffID = inferDiffID(cmd, inputnetworks);
		int consensusID = inferConsID(cmd, diffID);

		double cutoff = diffAlgo.default_weight_cutoff;
		if (cmd.hasOption(DiffanyOptions.cutoffShort))
		{
			cutoff = Double.parseDouble(cmd.getOptionValue(DiffanyOptions.cutoffShort));
		}

		/** THE ACTUAL ALGORITHM **/
		boolean cleanInput = true;
		Integer runID = p.addRunConfiguration(refNet, condNet, cleanInput, listener);
		Logger l = p.getLogger(runID);

		l.log("Calculating the pair-wise comparison between " + refNet.getName() + " and " + condNet.getName());

		// TODO v2.1: allow to change mode pairwise vs. differential
		diffAlgo.calculateOneDifferentialNetwork(p, runID, name, diffID, consensusID, cutoff, true, listener);

		// TODO v2.1: check number of differential networks generated
		RunOutput output = p.getOutput(runID);
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork diffNet = pair.getDifferentialNetwork();
		ConsensusNetwork consensusNet = pair.getConsensusNetwork();

		/** WRITE NETWORK OUTPUT **/
		boolean writeHeaders = true;

		File diffDir = getRequiredDir(cmd, DiffanyOptions.diffShort);
		NetworkIO.writeNetworkToDir(diffNet, diffDir, writeHeaders);
		l.log("Writing the differential network to " + diffDir);

		File consensusDir = getRequiredDir(cmd, DiffanyOptions.consensusShort);
		NetworkIO.writeNetworkToDir(consensusNet, consensusDir, writeHeaders);
		l.log("Writing the consensus network to " + consensusDir);

		l.log("Done !");

		/** WRITE LOG OUTPUT **/
		if (toLog)
		{
			for (LogEntry msg : l.getAllLogMessages())
			{
				System.out.println(msg);
			}
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
		if (cmd.hasOption(DiffanyOptions.consensusID))
		{
			consID = Integer.parseInt(cmd.getOptionValue(DiffanyOptions.consensusID));
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
