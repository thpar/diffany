package be.svlandeg.diffany.console;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.io.NetworkIO;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.OverlappingNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.DifferentialOutput;
import be.svlandeg.diffany.core.project.LogEntry;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.core.semantics.NodeMapper;

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
		
		// TODO v3.0: make ontologies adjustable
		Project p = new Project("Diffany-Analysis", new DefaultEdgeOntology(), new DefaultNodeMapper());
		NodeMapper nm = p.getNodeMapper();

		/** PARSE INPUT **/
		boolean toLog = false;
		if (cmd.hasOption(DiffanyOptions.logShort))
		{
			toLog = true;
		}

		String name = cmd.getOptionValue(DiffanyOptions.diffnameShort);

		double cutoff = diffAlgo.default_cutoff;
		if (cmd.hasOption(DiffanyOptions.cutoffShort))
		{
			cutoff = Double.parseDouble(cmd.getOptionValue(DiffanyOptions.cutoffShort));
		}
		
		File refDir = getRequiredDir(cmd, DiffanyOptions.refShort);
		ReferenceNetwork refNet = NetworkIO.readReferenceNetworkFromDir(refDir, nm);

		File condDir = getRequiredDir(cmd, DiffanyOptions.conShort);
		ConditionNetwork condNet = NetworkIO.readConditionNetworkFromDir(condDir, nm);

		/** THE ACTUAL ALGORITHM **/
		Integer rcID = p.addRunConfiguration(refNet, condNet);
		Logger l = p.getLogger(rcID);
		
		l.log("Calculating the pair-wise comparison between " + refNet.getName() + " and " + condNet.getName());
		
		// TODO v2.0: allow to change mode
		diffAlgo.calculateOneDifferentialNetwork(p, rcID, name, cutoff, true, true);
		
		// TODO v2.0: check number of differential networks generated
		DifferentialOutput output = p.getRunConfiguration(rcID).getDifferentialOutput();
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork diffNet = pair.getDifferentialNetwork();
		OverlappingNetwork overlapNet = pair.getOverlappingNetwork();

		/** WRITE NETWORK OUTPUT **/
		boolean allowVirtualEdges = true;
		
		File diffDir = getRequiredDir(cmd, DiffanyOptions.diffShort);
		NetworkIO.writeNetworkToDir(diffNet, nm, diffDir, allowVirtualEdges);
		l.log("Writing the differential network to " + diffDir);

		File overlapDir = getRequiredDir(cmd, DiffanyOptions.overlapShort);
		NetworkIO.writeNetworkToDir(overlapNet, nm, overlapDir, allowVirtualEdges);
		l.log("Writing the overlap network to " + overlapDir);
		
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
