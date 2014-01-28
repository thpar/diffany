package be.svlandeg.diffany.junit;

import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystem;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.OverlappingNetwork;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.examples.Bandyopadhyay2010;
import be.svlandeg.diffany.io.NetworkIO;

/** 
 * Class that automatically tests the IO functionality of Diffany:
 * writing/reading networks and projects to/from file.
 * 
 * TODO test reading
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
		String testLocation = System.getProperty("java.io.tmpdir") + "diffany" + File.separator;
		
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
		double cutoff = 0.0;
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, cutoff);
		
		// The reference network
		ReferenceNetwork rWriteNetwork = p.getReferenceNetwork();
		
		// The condition-dependent network (there is only 1 in that specific project)
		ConditionNetwork cWriteNetwork = p.getConditionNetworks().iterator().next();
		
		// There is exactly one differential network created
		Collection<DifferentialNetwork> dNetworks = p.getDifferentialNetworks();
		DifferentialNetwork dWriteNetwork = dNetworks.iterator().next();
				
		// The overlapping network (there should be only 1)
		OverlappingNetwork oWriteNetwork = dWriteNetwork.getOverlappingNetwork();
			
		// WRITING
		try
		{
			NetworkIO.writeReferenceNetworkToDir(rWriteNetwork, rDir);
			NetworkIO.writeConditionNetworkToDir(cWriteNetwork, cDir);
			NetworkIO.writeDifferentialNetworkToDir(dWriteNetwork, dDir);
			NetworkIO.writeOverlappingNetworkToDir(oWriteNetwork, oDir);
		}
		catch(IOException io)
		{
			fail(io.getMessage());
		}
		
		// READING
		try
		{
			ReferenceNetwork rReadNetwork = NetworkIO.readReferenceNetworkFromDir(rDir);
			ConditionNetwork cReadNetwork = NetworkIO.readConditionNetworkFromDir(cDir);
			Set<ConditionNetwork> cReadNetworks = new HashSet<ConditionNetwork>();
			cReadNetworks.add(cReadNetwork);
			
			DifferentialNetwork dReadNetwork = NetworkIO.readDifferentialNetworkFromDir(dDir, rReadNetwork, cReadNetworks);
			OverlappingNetwork oReadNetwork = NetworkIO.readOverlappingNetworkFromDir(oDir, rReadNetwork, cReadNetworks);
		}
		catch(IOException io)
		{
			fail(io.getMessage());
		}
		
	}
	
	
}
