package be.svlandeg.diffany.junit;

import static org.junit.Assert.*;

import java.util.Collection;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.OverlappingNetwork;
import be.svlandeg.diffany.core.project.DifferentialOutput;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.examples.MultipleConditionTest;

/**
 * Class that automatically tests the RunConfiguration settings and results stored in a Project entity.
 * This class only checks the type of networks generated, as other tests examine the content of those networks specifically.
 * 
 * @author Sofie Van Landeghem
 */
public class TestConfiguration
{

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results.
	 * 
	 * TODO check calc pairwise
	 */
	@Test
	public void testMultipleConditions()
	{
		CalculateDiff calc = new CalculateDiff();

		MultipleConditionTest ex = new MultipleConditionTest();
		double cutoff = 0.0;
		Project p = ex.getTestProject();
		assertEquals(0, p.getAllRunConfigurations().size());
		
		// CONFIG 1: 1-all differential configuration - before calculation
		int ID1 = ex.getTestDiffConfiguration(p);
		assertEquals(1, p.getAllRunConfigurations().size());
		Collection<DifferentialOutput> dOutputs1 = p.getRunConfiguration(ID1).getDifferentialOutputs();
		assertEquals(0, dOutputs1.size());
		
		// CONFIG 1: 1-all differential configuration : differential+overlap calculation
		calc.calculateOneDifferentialNetwork(p, ID1, cutoff, true, true);
		dOutputs1 = p.getRunConfiguration(ID1).getDifferentialOutputs();
		assertEquals(1, dOutputs1.size());

		DifferentialOutput output1 = dOutputs1.iterator().next();
		OutputNetworkPair pair1 = output1.getOutputAsPair();

		DifferentialNetwork dNetwork1 = pair1.getDifferentialNetwork();
		assertEquals(8, dNetwork1.getEdges().size());

		OverlappingNetwork oNetwork1 = pair1.getOverlappingNetwork();
		assertEquals(3, oNetwork1.getEdges().size());

		
		// CONFIG 2: (the exact same) 1-all differential+overlap calculation
		int ID2 = ex.getTestDiffConfiguration(p);
		calc.calculateOneDifferentialNetwork(p, ID2, cutoff, true, true);
		assertEquals(2, p.getAllRunConfigurations().size());

		Collection<DifferentialOutput> dOutputs2 = p.getRunConfiguration(ID2).getDifferentialOutputs();
		assertEquals(1, dOutputs2.size());

		DifferentialOutput output2 = dOutputs2.iterator().next();
		OutputNetworkPair pair2 = output2.getOutputAsPair();

		DifferentialNetwork dNetwork2 = pair2.getDifferentialNetwork();
		assertEquals(8, dNetwork2.getEdges().size());

		OverlappingNetwork oNetwork2 = pair2.getOverlappingNetwork();
		assertEquals(3, oNetwork2.getEdges().size());
		
		
		// CONFIG 3: (only) overlap configuration - assess correct failure of differential calculation
		int ID3 = ex.getTestOverlapConfiguration(p);
		assertEquals(3, p.getAllRunConfigurations().size());
		
		boolean exceptionThrown1 = false;
		try
		{
			calc.calculateOneDifferentialNetwork(p, ID3, cutoff, true, true);
		}
		catch(IllegalArgumentException e)
		{
			exceptionThrown1 = true;
		}
		assertTrue(exceptionThrown1);
		Collection<DifferentialOutput> dOutputs3 = p.getRunConfiguration(ID3).getDifferentialOutputs();
		assertEquals(0, dOutputs3.size());
		
		// CONFIG 3: (only) overlap calculation - assess correct failure of differential network request
		calc.calculateOneDifferentialNetwork(p, ID3, cutoff, false, true);
		dOutputs3 = p.getRunConfiguration(ID3).getDifferentialOutputs();
		assertEquals(1, dOutputs3.size());
		
		DifferentialOutput output3 = dOutputs3.iterator().next();
		
		boolean exceptionThrown2 = false;
		try
		{
			@SuppressWarnings("unused")
            OutputNetworkPair pair3 = output3.getOutputAsPair();
		}
		catch(IllegalArgumentException e)
		{
			exceptionThrown2 = true;
		}
		assertTrue(exceptionThrown2);

		//DifferentialNetwork dNetwork3 = pair3.getDifferentialNetwork();
		//assertEquals(8, dNetwork3.getEdges().size());

		OverlappingNetwork oNetwork3 = output3.getOverlappingNetwork();
		assertEquals(3, oNetwork3.getEdges().size());
	}

}
