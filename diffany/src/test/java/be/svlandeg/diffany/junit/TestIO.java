package be.svlandeg.diffany.junit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.io.NetworkIO;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.OverlappingNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunConfiguration;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.examples.Bandyopadhyay2010;

/** 
 * Class that automatically tests the IO functionality of Diffany:
 * writing/reading networks and projects to/from file.
 * 
 * @author Sofie Van Landeghem
 */
public class TestIO
{
	
	/**
	 * JUNIT Test: write and read all types of Networks.
	 * The networks used in this tests are tested for completeness in {@link TestExamples#testBandyopadhyay}
	 */
	@Test
	public void testNetworkIO()
	{
		// System-dependent tmp dir. E.g. windows 7: C:\Users\YourUserName\AppData\Local\Temp
		String testLocation = System.getProperty("java.io.tmpdir") + File.separator + "diffany" + File.separator;
		
		File rDir = new File(testLocation + "reference/");
		File cDir = new File(testLocation + "condition/");
		File dDir = new File(testLocation + "differential/");
		File oDir = new File(testLocation + "overlap/");
		
		rDir.mkdirs();
		cDir.mkdirs();
		dDir.mkdirs();
		oDir.mkdirs();
		
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		Project p = ex.getProjectFigure1C();
		int ID = ex.getTestConfiguration1C(p);
		double cutoff = 0.0;
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, ID, cutoff);
		
		RunConfiguration rc = p.getRunConfiguration(ID);
		
		// The reference network
		ReferenceNetwork rWriteNetwork = rc.getReferenceNetwork();
		
		// The condition-dependent network (there is only 1 in that specific project)
		ConditionNetwork cWriteNetwork = rc.getConditionNetworks().iterator().next();
		
		// There is exactly one differential network created
		Collection<DifferentialNetwork> dNetworks = rc.getDifferentialNetworks();
		DifferentialNetwork dWriteNetwork = dNetworks.iterator().next();
				
		// The overlapping network (there should be only 1)
		OverlappingNetwork oWriteNetwork = dWriteNetwork.getOverlappingNetwork();
		
		NodeMapper nm = new DefaultNodeMapper();
			
		// WRITING
		try
		{
			NetworkIO.writeReferenceNetworkToDir(rWriteNetwork, nm, rDir);
			NetworkIO.writeConditionNetworkToDir(cWriteNetwork, nm, cDir);
			NetworkIO.writeDifferentialNetworkToDir(dWriteNetwork, nm, dDir);
			NetworkIO.writeOverlappingNetworkToDir(oWriteNetwork, nm, oDir);
		}
		catch(IOException io)
		{
			fail(io.getMessage());
		}
		
		// READING
		try
		{
			ReferenceNetwork rReadNetwork = NetworkIO.readReferenceNetworkFromDir(rDir, nm);
			assertEquals(3, rReadNetwork.getEdges().size());
			assertEquals(5, rReadNetwork.getNodes().size());
			
			ConditionNetwork cReadNetwork = NetworkIO.readConditionNetworkFromDir(cDir, nm);
			assertEquals(3, cReadNetwork.getEdges().size());
			assertEquals(4, cReadNetwork.getNodes().size());
			
			Set<ConditionNetwork> cReadNetworks = new HashSet<ConditionNetwork>();
			cReadNetworks.add(cReadNetwork);
			
			DifferentialNetwork dReadNetwork = NetworkIO.readDifferentialNetworkFromDir(dDir, nm, rReadNetwork, cReadNetworks);
			assertEquals(3, dReadNetwork.getEdges().size());
			assertEquals(4, dReadNetwork.getNodes().size());
			
			OverlappingNetwork oReadNetwork = NetworkIO.readOverlappingNetworkFromDir(oDir, nm, rReadNetwork, cReadNetworks);
			assertEquals(2, oReadNetwork.getEdges().size());
			assertEquals(3, oReadNetwork.getNodes().size());
			
			
		}
		catch(IOException io)
		{
			fail(io.getMessage());
		}
		
	}
	
	
}
