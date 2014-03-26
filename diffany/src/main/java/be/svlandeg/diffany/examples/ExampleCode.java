package be.svlandeg.diffany.examples;

import java.io.File;
import java.io.IOException;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.Logger;
import be.svlandeg.diffany.core.Project;
import be.svlandeg.diffany.core.RunConfiguration;
import be.svlandeg.diffany.io.NetworkIO;
import be.svlandeg.diffany.networks.ConditionNetwork;
import be.svlandeg.diffany.networks.DifferentialNetwork;
import be.svlandeg.diffany.networks.OverlappingNetwork;
import be.svlandeg.diffany.networks.ReferenceNetwork;
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

		File condDir = new File(condLocation);
		ConditionNetwork condNet = NetworkIO.readConditionNetworkFromDir(condDir, nm);

		/** DEFINE THE RUN PARAMETERS **/
		double cutoff = 0.0;
		int rcID = p.addRunConfiguration(refNet, condNet);
		
		/** THE ACTUAL ALGORITHM **/
		CalculateDiff diffAlgo = new CalculateDiff();
		diffAlgo.calculateOneDifferentialNetwork(p, rcID, cutoff);

		// In this case, there will be exactly one DifferentialNetwork
		RunConfiguration rc = p.getRunConfiguration(rcID);
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
