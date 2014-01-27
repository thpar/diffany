package be.svlandeg.diffany.junit;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.Collection;

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
 * @author Sofie Van Landeghem
 */
public class TestIO
{
	
	private static String testLocation = "C:/temp/diffany/";
	
	/**
	 * JUNIT Test: write and read a {@link ReferenceNetwork}.
	 * The networks used in this tests are tested for completeness in {@link TestExamples#testBandyopadhyay}
	 */
	@Test
	public void testReferenceNetworkIO()
	{
		File testDir = new File(testLocation + "reference/");
		testDir.mkdirs();
		
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		Project p = ex.getProjectFigure1C();
		
		// The reference network
		ReferenceNetwork rNetwork = p.getReferenceNetwork();
		boolean written = true;
		try
		{
			NetworkIO.writeNetworkToDir(rNetwork, testDir);
		}
		catch(IOException io)
		{
			written = false;
		}
		assertTrue(written);
	}
	
	/**
	 * JUNIT Test: write and read a {@link ConditionNetwork}.
	 * The networks used in this tests are tested for completeness in {@link TestExamples#testBandyopadhyay}
	 */
	@Test
	public void testConditionNetworkIO()
	{
		File testDir = new File(testLocation + "condition");
		testDir.mkdirs();
		
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		double cutoff = 0.0;
		Project p = ex.getProjectFigure1C();
		
		// The condition-dependent network (there is only 1 in that specific project)
		ConditionNetwork cNetwork = p.getConditionNetworks().iterator().next();
		boolean written = true;
		try
		{
			NetworkIO.writeNetworkToDir(cNetwork, testDir);
		}
		catch(IOException io)
		{
			written = false;
		}
		assertTrue(written);
	}

	
	/**
	 * JUNIT Test: write and read a {@link DifferentialNetwork}.
	 * The networks used in this tests are tested for completeness in {@link TestExamples#testBandyopadhyay}
	 */
	@Test
	public void testDifferentialNetworkIO()
	{
		File testDir = new File(testLocation + "differential/");
		testDir.mkdirs();
		
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		double cutoff = 0.0;
		Project p = ex.getProjectFigure1C();
		
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, cutoff);
		
		// There is exactly one differential network created
		Collection<DifferentialNetwork> dNetworks = p.getDifferentialNetworks();
		DifferentialNetwork dNetwork = dNetworks.iterator().next();
		
		boolean written = true;
		try
		{
			NetworkIO.writeNetworkToDir(dNetwork, testDir);
		}
		catch(IOException io)
		{
			written = false;
		}
		assertTrue(written);
	}
	
	/**
	 * JUNIT Test: write and read a {@link OverlappingNetwork}.
	 * The networks used in this tests are tested for completeness in {@link TestExamples#testBandyopadhyay}
	 */
	@Test
	public void testOverlappingNetworkIO()
	{
		File testDir = new File(testLocation + "overlap/");
		testDir.mkdirs();
		
		Bandyopadhyay2010 ex = new Bandyopadhyay2010();
		double cutoff = 0.0;
		Project p = ex.getProjectFigure1C();
		
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(p, cutoff);
		
		// There is exactly one differential network created
		Collection<DifferentialNetwork> dNetworks = p.getDifferentialNetworks();
		DifferentialNetwork dNetwork = dNetworks.iterator().next();
		
		// The overlapping network (there should be only 1)
		OverlappingNetwork sNetwork = dNetwork.getOverlappingNetwork();
		
		boolean written = true;
		try
		{
			NetworkIO.writeNetworkToDir(sNetwork, testDir);
		}
		catch(IOException io)
		{
			written = false;
		}
		assertTrue(written);
	}
	
	
}
