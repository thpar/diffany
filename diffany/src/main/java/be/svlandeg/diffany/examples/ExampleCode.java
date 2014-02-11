package be.svlandeg.diffany.examples;

import java.io.File;
import java.io.IOException;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Logger;
import be.svlandeg.diffany.concepts.OverlappingNetwork;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.concepts.RunConfiguration;
import be.svlandeg.diffany.io.NetworkIO;
import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.semantics.NodeMapper;

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
		EdgeOntology eo = new DefaultEdgeOntology();
		NodeMapper nm = new DefaultNodeMapper();
		Project p = new Project("testProject", eo, nm);

		/** READ THE INPUT NETWORKS **/
		File refDir = new File(refLocation);
		ReferenceNetwork refNet = NetworkIO.readReferenceNetworkFromDir(refDir, nm);
		p.registerSourceNetwork(refNet);

		File condDir = new File(condLocation);
		ConditionNetwork condNet = NetworkIO.readConditionNetworkFromDir(condDir, nm);
		p.registerSourceNetwork(condNet);

		/** DEFINE THE RUN PARAMETERS **/
		double cutoff = 0.0;
		RunConfiguration rc = new RunConfiguration(refNet, condNet);
		int rcID = p.addRunConfiguration(rc);
		
		/** THE ACTUAL ALGORITHM **/
		CalculateDiff diffAlgo = new CalculateDiff();
		diffAlgo.calculateOneDifferentialNetwork(p, rcID, cutoff);

		// In this case, there will be exactly one DifferentialNetwork
		DifferentialNetwork diffNet = rc.getDifferentialNetworks().iterator().next();
		OverlappingNetwork overlapNet = diffNet.getOverlappingNetwork();

		/** WRITE NETWORK OUTPUT **/
		File diffDir = new File(diffLocation);
		NetworkIO.writeDifferentialNetworkToDir(diffNet, nm, diffDir);

		File overlapDir = new File(overlapLocation);
		NetworkIO.writeOverlappingNetworkToDir(overlapNet, nm, overlapDir);

		/** WRITE LOG OUTPUT **/
		Logger logger = p.getLogger(rcID);
		for (String msg : logger.getAllLogMessages())
		{
			System.out.println(msg);
		}
	}

}
