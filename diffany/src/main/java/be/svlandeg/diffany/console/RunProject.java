package be.svlandeg.diffany.console;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Logger;
import be.svlandeg.diffany.concepts.OverlappingNetwork;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.io.NetworkIO;
import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

/**
 * This class can run the Diffany algorithms from a {@link org.apache.commons.cli.CommandLine} object.
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

		Logger logger = new Logger();

		// TODO v2.0: adjustable
		EdgeOntology eo = new DefaultEdgeOntology();

		// TODO v2.0: adjustable
		NodeMapper nm = new DefaultNodeMapper();

		File refDir = getRequiredDir(cmd, DiffanyOptions.refShort);
		ReferenceNetwork refNet = NetworkIO.readReferenceNetworkFromDir(refDir, nm);

		File condDir = getRequiredDir(cmd, DiffanyOptions.conShort);
		ConditionNetwork condNet = NetworkIO.readConditionNetworkFromDir(condDir, nm);

		/** THE ACTUAL ALGORITHM **/
		logger.log("Calculating the pair-wise comparison between " + refNet.getName() + " and " + condNet.getName());
		DifferentialNetwork diffNet = diffAlgo.calculateDiffNetwork(refNet, condNet, eo, nm, name, cutoff, logger);
		OverlappingNetwork overlapNet = diffNet.getOverlappingNetwork();

		/** WRITE NETWORK OUTPUT **/
		File diffDir = getRequiredDir(cmd, DiffanyOptions.diffShort);
		NetworkIO.writeDifferentialNetworkToDir(diffNet, nm, diffDir);
		logger.log("Writing the differential network to " + diffDir);

		File overlapDir = getRequiredDir(cmd, DiffanyOptions.overlapShort);
		NetworkIO.writeOverlappingNetworkToDir(overlapNet, nm, overlapDir);
		logger.log("Writing the overlap network to " + overlapDir);
		
		logger.log("Done !");

		/** WRITE LOG OUTPUT **/
		if (toLog)
		{
			for (String msg : logger.getAllLogMessages())
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
