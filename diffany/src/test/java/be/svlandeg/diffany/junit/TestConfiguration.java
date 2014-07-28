package be.svlandeg.diffany.junit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.OverlappingNetwork;
import be.svlandeg.diffany.core.project.RunOutput;
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

	private static double cutoff = 0.0;

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results.
	 */
	@Test
	public void testMultipleConditions()
	{
		CalculateDiff calc = new CalculateDiff();

		MultipleConditionTest ex = new MultipleConditionTest();
		Project p = ex.getTestProject();

		int calls = 0;

		// before calculation
		testIni(p, ex, calls);
		calls++;

		// 1-all differential configuration : differential+overlap calculation
		testDifferentialOneMulti(p, calc, ex, calls, 10, 11);
		calls++;

		// (the exact same) 1-all differential+overlap calculation
		testDifferentialOneMulti(p, calc, ex, calls, 12, 13);
		calls++;

		// 1-all differential (only) calculation
		testDifferentialOneMulti(p, calc, ex, calls, 14, -1);
		calls++;

		// 1-all overlap (only) calculation - with differential configuration
		testDifferentialOneMulti(p, calc, ex, calls, -1, 15);
		calls++;

		// 1-all calculation which doesn't actually do anything (no diff, no overlap)
		testDifferentialOneMulti(p, calc, ex, calls, -1, -1);
		calls++;

		// 1-all overlap configuration - assess correct failure of differential calculation
		testOverlapOneMulti(p, calc, ex, calls, 20);
		calls++;

		// pairwise differential configuration : differential+overlap calculation
		testDifferentialPairwise(p, calc, ex, calls, true, true, 30);
		calls++;

		// pairwise differential (only) calculation
		testDifferentialPairwise(p, calc, ex, calls, true, false, 40);
		calls++;

		// pairwise overlap (only) calculation - with differential configuration
		testDifferentialPairwise(p, calc, ex, calls, false, false, 50);
		calls++;

		// pairwise calculation which doesn't actually do anything (no diff, no overlap)
		testDifferentialPairwise(p, calc, ex, calls, false, false, 60);
		calls++;
	}

	/**
	 * Test an initial empty project: should contain no differential output.
	 */
	private void testIni(Project p, MultipleConditionTest ex, int calls)
	{
		assertEquals(0, p.getAllRunIDs().size());

		int ID = ex.getTestDiffConfiguration(p);
		assertEquals(calls + 1, p.getAllRunIDs().size());
		
		RunOutput dOutput = p.getOutput(ID);
		assertNrDiffNetworks(dOutput, 0);
	}

	/**
	 * Test a 1-multi differential configuration, with differential and/or overlap calculation.
	 */
	private void testDifferentialOneMulti(Project p, CalculateDiff calc, MultipleConditionTest ex, int calls, int diffID, int overlapID)
	{
		int ID = ex.getTestDiffConfiguration(p);
		assertEquals(calls + 1, p.getAllRunIDs().size());

		// it should not make a difference how many times this method is called!
		calc.calculateOneDifferentialNetwork(p, ID, cutoff, diffID, overlapID, true);
		calc.calculateOneDifferentialNetwork(p, ID, cutoff, diffID, overlapID, true);
		RunOutput output = p.getOutput(ID);
		
		boolean diff = diffID >= 0;
		boolean overlap = overlapID >= 0;
		
		if (diff && overlap)
		{
			assertNrPairs(output, 1);
			OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
			
			DifferentialNetwork dNetwork = pair.getDifferentialNetwork();
			assertEquals(8, dNetwork.getEdges().size());

			OverlappingNetwork oNetwork = pair.getOverlappingNetwork();
			assertEquals(3, oNetwork.getEdges().size());
		}
		else
		{
			assertNrPairs(output, 0);
		}
		if (diff)
		{
			assertNrDiffNetworks(output, 1);
			DifferentialNetwork dNetwork = output.getDifferentialNetworks().iterator().next();
			assertEquals(8, dNetwork.getEdges().size());
		}
		else
		{
			assertNrDiffNetworks(output, 0);
		}
		if (overlap)
		{
			assertNrOverlapNetworks(output, 1);
			assertEquals(1, output.getOverlappingNetworks().size());
			OverlappingNetwork oNetwork = output.getOverlappingNetworks().iterator().next();
			assertEquals(3, oNetwork.getEdges().size());
		}
		else
		{
			assertNrOverlapNetworks(output, 0);
		}
	}

	/**
	 * Test a 1-multi overlap configuration, should not contain differential networks.
	 */
	private void testOverlapOneMulti(Project p, CalculateDiff calc, MultipleConditionTest ex, int calls, int overlapID)
	{
		int ID = ex.getTestOverlapConfiguration(p);
		assertEquals(calls + 1, p.getAllRunIDs().size());

		assertException(p, calc, ID, overlapID++, overlapID++);
		assertException(p, calc, ID, overlapID++, -1);

		RunOutput dOutput = p.getOutput(ID);
		assertNrPairs(dOutput, 0);
		assertNrDiffNetworks(dOutput, 0);
		assertNrOverlapNetworks(dOutput, 0);

		calc.calculateOneDifferentialNetwork(p, ID, cutoff, -1, overlapID, true);
		dOutput = p.getOutput(ID);

		assertNrPairs(dOutput, 0);
		assertNrDiffNetworks(dOutput, 0);
		assertNrOverlapNetworks(dOutput, 1);

		OverlappingNetwork oNetwork = dOutput.getOverlappingNetworks().iterator().next();
		assertEquals(3, oNetwork.getEdges().size());
	}

	/**
	 * Test a pairwise differential configuration, with differential and/or overlap calculation.
	 */
	private void testDifferentialPairwise(Project p, CalculateDiff calc, MultipleConditionTest ex, int calls, boolean diff, boolean overlap, int firstID)
	{
		int ID = ex.getTestDiffConfiguration(p);
		assertEquals(calls + 1, p.getAllRunIDs().size());

		calc.calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, diff, overlap, firstID, true);
		RunOutput output = p.getOutput(ID);

		// First, test the number of result networks, depending on the settings.
		if (diff)
		{
			// ref vs. 2 conditions
			assertNrDiffNetworks(output, 2);
			if (overlap)
			{
				assertNrOverlapNetworks(output, 2);
			}
			else
			{
				assertNrOverlapNetworks(output, 0);
			}
		}
		
		else if (overlap)
		{
			// 3 pairwise comparisons
			assertNrPairs(output, 0);
			assertNrDiffNetworks(output, 0);
			assertNrOverlapNetworks(output, 3);
		}
		else
		{
			// no comparisons
			assertNrPairs(output, 0);
			assertNrDiffNetworks(output, 0);
			assertNrOverlapNetworks(output, 0);
		}

		// Now, test the number of edges per result network, depending on the settings.
		Map<String, Integer> diffcounts = new HashMap<String, Integer>();
		diffcounts.put("diff_Draughty", 12);
		diffcounts.put("diff_Salty", 11);

		Map<String, Integer> overlapcounts = new HashMap<String, Integer>();
		overlapcounts.put("overlap_Draughty_Reference", 3);
		overlapcounts.put("overlap_Reference_Salty", 5);
		overlapcounts.put("overlap_Draughty_Salty", 5);	

		if (diff && overlap)
		{
			for (OutputNetworkPair pair : output.getOutputAsPairs())
			{
				DifferentialNetwork dNetwork = pair.getDifferentialNetwork();
				int countD = diffcounts.get(dNetwork.getName());
				assertEquals(countD, dNetwork.getEdges().size());

				OverlappingNetwork oNetwork = pair.getOverlappingNetwork();
				int countO = overlapcounts.get(oNetwork.getName());
				assertEquals(countO, oNetwork.getEdges().size());
			}
		}
		if (diff)
		{
			for (DifferentialNetwork dNetwork : output.getDifferentialNetworks())
			{
				int countD = diffcounts.get(dNetwork.getName());
				assertEquals(countD, dNetwork.getEdges().size());
			}
		}
		if (overlap)
		{
			for (OverlappingNetwork oNetwork : output.getOverlappingNetworks())
			{
				int countO = overlapcounts.get(oNetwork.getName());
				assertEquals(countO, oNetwork.getEdges().size());
			}
		}
	}

	/**
	 * Private method that asserts an error will result from calling the oneMulti calculation with a specific runconfiguration ID.
	 */
	private void assertException(Project p, CalculateDiff calc, int ID, int diffID, int overlapID)
	{
		boolean exceptionThrown = false;
		try
		{
			calc.calculateOneDifferentialNetwork(p, ID, cutoff, diffID, overlapID, true);
		}
		catch (IllegalArgumentException e)
		{
			exceptionThrown = true;
		}
		assertTrue(exceptionThrown);
	}

	/**
	 * Private method that asserts the number of differential output pairs in the output result (may be 0).
	 */
	private void assertNrPairs(RunOutput output, int number)
	{
		int pairs = output.getOutputAsPairs().size();
		assertEquals(number, pairs);
	}

	/**
	 * Private method that asserts the number of differential networks in the output result (may be 0).
	 */
	private void assertNrDiffNetworks(RunOutput output, int number)
	{
		int diffs = output.getDifferentialNetworks().size();
		assertEquals(number, diffs);
	}

	/**
	 * Private method that asserts the number of overlap networks in the output result (may be 0).
	 */
	private void assertNrOverlapNetworks(RunOutput output, int number)
	{
		int overlaps = output.getOverlappingNetworks().size();
		assertEquals(number, overlaps);
	}

}
