package be.svlandeg.diffany.examples;

import java.io.File;
import java.io.IOException;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.io.NetworkIO;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.core.project.LogEntry;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

/** 
 * This class provides example code to work with the Diffany library.
 * 
 * @author Sofie Van Landeghem
 */
public class ExampleCode
{
	
	public void example(String refLocation, String condLocation, String diffLocation, String consensusLocation) throws IOException
	{
		/** DEFINE THE ONTOLOGIES AND THE PROJECT **/
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project("testProject", eo);

		/** READ THE INPUT NETWORKS **/
		boolean skipHeader = true;
		File refDir = new File(refLocation);
		ReferenceNetwork refNet = NetworkIO.readReferenceNetworkFromDir(refDir, skipHeader);

		File condDir = new File(condLocation);
		ConditionNetwork condNet = NetworkIO.readConditionNetworkFromDir(condDir, skipHeader);

		/** DEFINE THE RUN PARAMETERS **/
		double cutoff = 0.0;
		boolean cleanInput = true;
		int runID = p.addRunConfiguration(refNet, condNet, cleanInput, null);
		
		/** THE ACTUAL ALGORITHM **/
		CalculateDiff diffAlgo = new CalculateDiff();
		diffAlgo.calculateOneDifferentialNetwork(p, runID, cutoff, null, null, 342, 666, true, null);

		// In this case, there will be exactly one DifferentialNetwork
		RunOutput output = p.getOutput(runID);
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork diffNet = pair.getDifferentialNetwork();
		ConsensusNetwork consensusNet = pair.getConsensusNetwork();

		/** WRITE NETWORK OUTPUT **/
		boolean writeHeaders = true;
		
		File diffDir = new File(diffLocation);
		NetworkIO.writeNetworkToDir(diffNet, diffDir, writeHeaders);

		File consensusDir = new File(consensusLocation);
		NetworkIO.writeNetworkToDir(consensusNet, consensusDir, writeHeaders);

		/** WRITE LOG OUTPUT **/
		Logger logger = p.getLogger(runID);
		for (LogEntry msg : logger.getAllLogMessages())
		{
			System.out.println(msg);
		}
	}

}
