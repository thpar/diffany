package be.svlandeg.diffany.examples;

import java.io.File;
import java.io.IOException;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.io.NetworkIO;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.OverlappingNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.DifferentialOutput;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunConfiguration;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

/** 
 * This class provides example code to work with the Diffany library.
 * 
 * @author Sofie Van Landeghem
 */
public class ExampleCode
{
	
	public void example(String refLocation, String condLocation, String diffLocation, String overlapLocation) throws IOException
	{
		/** DEFINE THE ONTOLOGIES AND THE PROJECT **/
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		NodeMapper nm = new DefaultNodeMapper();
		Project p = new Project("testProject", eo, nm);

		/** READ THE INPUT NETWORKS **/
		File refDir = new File(refLocation);
		ReferenceNetwork refNet = NetworkIO.readReferenceNetworkFromDir(refDir, nm);

		File condDir = new File(condLocation);
		ConditionNetwork condNet = NetworkIO.readConditionNetworkFromDir(condDir, nm);

		/** DEFINE THE RUN PARAMETERS **/
		double cutoff = 0.0;
		int rcID = p.addRunConfiguration(refNet, condNet);
		
		/** THE ACTUAL ALGORITHM **/
		CalculateDiff diffAlgo = new CalculateDiff();
		diffAlgo.calculateOneDifferentialNetwork(p, rcID, cutoff, true, true);

		// In this case, there will be exactly one DifferentialNetwork
		RunConfiguration rc = p.getRunConfiguration(rcID);
		DifferentialOutput output = rc.getDifferentialOutput();
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
		DifferentialNetwork diffNet = pair.getDifferentialNetwork();
		OverlappingNetwork overlapNet = pair.getOverlappingNetwork();

		/** WRITE NETWORK OUTPUT **/
		File diffDir = new File(diffLocation);
		NetworkIO.writeNetworkToDir(diffNet, nm, diffDir);

		File overlapDir = new File(overlapLocation);
		NetworkIO.writeNetworkToDir(overlapNet, nm, overlapDir);

		/** WRITE LOG OUTPUT **/
		Logger logger = p.getLogger(rcID);
		for (String msg : logger.getAllLogMessages())
		{
			System.out.println(msg);
		}
	}

}
