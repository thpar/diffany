package be.svlandeg.diffany.junit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

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

	private static double cutoff = 0.0;

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results.
	 * 
	 * TODO check calc pairwise
	 */
	//@Test
	public void testMultipleConditions()
	{
		CalculateDiff calc = new CalculateDiff();

		MultipleConditionTest ex = new MultipleConditionTest();
		Project p = ex.getTestProject();

		int calls = 0;

		// before calculation
		testIni(p, calc, ex, calls);
		calls++;

		// 1-all differential configuration : differential+overlap calculation
		testDifferentialOneMulti(p, calc, ex, calls, true, true);
		calls++;

		// (the exact same) 1-all differential+overlap calculation
		testDifferentialOneMulti(p, calc, ex, calls, true, true);
		calls++;

		// 1-all differential (only) calculation
		testDifferentialOneMulti(p, calc, ex, calls, true, false);
		calls++;

		// 1-all overlap (only) calculation - with differential configuration
		testDifferentialOneMulti(p, calc, ex, calls, false, true);
		calls++;

		// 1-all calculation which doesn't actually do anything (no diff, no overlap)
		testDifferentialOneMulti(p, calc, ex, calls, false, false);
		calls++;

		// 1-all overlap configuration - assess correct failure of differential calculation
		testOverlapOneMulti(p, calc, ex, calls);
		calls++;

		// pairwise differential configuration : differential+overlap calculation
		testDifferentialPairwise(p, calc, ex, calls, true, true);
		calls++;
		
		// pairwise differential (only) calculation
		testDifferentialPairwise(p, calc, ex, calls, true, false);
		calls++;
				
		// pairwise overlap (only) calculation - with differential configuration
		testDifferentialPairwise(p, calc, ex, calls, false, true);
		calls++;
		
		// pairwise calculation which doesn't actually do anything (no diff, no overlap)
		testDifferentialPairwise(p, calc, ex, calls, false, false);
		calls++;
	}

	/**
	 * Test an initial empty project 
	 */
	private void testIni(Project p, CalculateDiff calc, MultipleConditionTest ex, int calls)
	{
		assertEquals(0, p.getAllRunConfigurations().size());

		int ID = ex.getTestDiffConfiguration(p);
		assertEquals(calls + 1, p.getAllRunConfigurations().size());
		Collection<DifferentialOutput> dOutputs = p.getRunConfiguration(ID).getDifferentialOutputs();
		assertEquals(0, dOutputs.size());
	}

	/**
	 * Test an 1-multi differential configuration
	 */
	private void testDifferentialOneMulti(Project p, CalculateDiff calc, MultipleConditionTest ex, int calls, boolean diff, boolean overlap)
	{
		int ID = ex.getTestDiffConfiguration(p);
		assertEquals(calls + 1, p.getAllRunConfigurations().size());

		// it should not make a difference how many times this method is called!
		calc.calculateOneDifferentialNetwork(p, ID, cutoff, diff, overlap);
		calc.calculateOneDifferentialNetwork(p, ID, cutoff, diff, overlap);
		Collection<DifferentialOutput> dOutputs = p.getRunConfiguration(ID).getDifferentialOutputs();
		dOutputs = p.getRunConfiguration(ID).getDifferentialOutputs();
		assertEquals(1, dOutputs.size());

		DifferentialOutput output = dOutputs.iterator().next();

		if (diff && overlap)
		{
			OutputNetworkPair pair = output.getOutputAsPair();
			DifferentialNetwork dNetwork = pair.getDifferentialNetwork();
			assertEquals(8, dNetwork.getEdges().size());

			OverlappingNetwork oNetwork = pair.getOverlappingNetwork();
			assertEquals(3, oNetwork.getEdges().size());
		}
		else
		{
			assertNoPair(output);
		}
		if (diff)
		{
			DifferentialNetwork dNetwork = output.getDifferentialNetwork();
			assertEquals(8, dNetwork.getEdges().size());
		}
		else
		{
			assertNoDiffNetwork(output);
		}
		if (overlap)
		{
			OverlappingNetwork oNetwork5 = output.getOverlappingNetwork();
			assertEquals(3, oNetwork5.getEdges().size());
		}
		else
		{
			assertNoOverlapNetwork(output);
		}
	}

	/**
	 * Test an 1-multi overlap configuration
	 */
	private void testOverlapOneMulti(Project p, CalculateDiff calc, MultipleConditionTest ex, int calls)
	{
		int ID = ex.getTestOverlapConfiguration(p);
		assertEquals(calls + 1, p.getAllRunConfigurations().size());

		assertException(p, calc, ID, true, true);
		assertException(p, calc, ID, true, false);

		Collection<DifferentialOutput> dOutputs = p.getRunConfiguration(ID).getDifferentialOutputs();
		assertEquals(0, dOutputs.size());

		calc.calculateOneDifferentialNetwork(p, ID, cutoff, false, true);
		dOutputs = p.getRunConfiguration(ID).getDifferentialOutputs();
		assertEquals(1, dOutputs.size());

		DifferentialOutput output = dOutputs.iterator().next();

		assertNoPair(output);
		assertNoDiffNetwork(output);

		OverlappingNetwork oNetwork = output.getOverlappingNetwork();
		assertEquals(3, oNetwork.getEdges().size());
	}

	/**
	 * Test an 1-multi differential configuration
	 */
	private void testDifferentialPairwise(Project p, CalculateDiff calc, MultipleConditionTest ex, int calls, boolean diff, boolean overlap)
	{
		int ID = ex.getTestDiffConfiguration(p);
		assertEquals(calls + 1, p.getAllRunConfigurations().size());

		calc.calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, diff, overlap);
		Collection<DifferentialOutput> dOutputs = p.getRunConfiguration(ID).getDifferentialOutputs();
		dOutputs = p.getRunConfiguration(ID).getDifferentialOutputs();
		
		if (diff)
		{
			// ref vs. 2 conditions
			assertEquals(2, dOutputs.size());
		}
		else if (overlap)
		{
			// 3 pairwise comparisons
			assertEquals(3, dOutputs.size());
		}
		else
		{
			// no comparisons
			assertEquals(0, dOutputs.size());
		}

		for (DifferentialOutput output : dOutputs)
		{
			Map<String, Integer> diffcounts = new HashMap<String, Integer>();
			diffcounts.put("diff_Draughty", 12);
			diffcounts.put("diff_Salty", 11);
			
			Map<String, Integer> overlapcounts = new HashMap<String, Integer>();
			overlapcounts.put("overlap_Reference_Draughty", 3);
			overlapcounts.put("overlap_Draughty_Reference", 3);
			overlapcounts.put("overlap_Reference_Salty", 5);
			overlapcounts.put("overlap_Salty_Reference", 5);
			overlapcounts.put("overlap_Draughty_Salty", 5);
			overlapcounts.put("overlap_Salty_Draughty", 5);
			
			
			if (diff && overlap)
			{
				OutputNetworkPair pair = output.getOutputAsPair();
				DifferentialNetwork dNetwork = pair.getDifferentialNetwork();
				int countD = diffcounts.get(dNetwork.getName());
				assertEquals(countD, dNetwork.getEdges().size());

				OverlappingNetwork oNetwork = pair.getOverlappingNetwork();
				int countO = overlapcounts.get(oNetwork.getName());
				System.out.println(oNetwork.getName() + " - " + countO);
				System.out.println("found " + oNetwork.getEdges().size());
				assertEquals(countO, oNetwork.getEdges().size());
			}
			else
			{
				assertNoPair(output);
			}
			if (diff)
			{
				DifferentialNetwork dNetwork = output.getDifferentialNetwork();
				int countD = diffcounts.get(dNetwork.getName());
				assertEquals(countD, dNetwork.getEdges().size());
			}
			else
			{
				assertNoDiffNetwork(output);
			}
			if (overlap)
			{
				OverlappingNetwork oNetwork = output.getOverlappingNetwork();
				int countO = overlapcounts.get(oNetwork.getName());
				assertEquals(countO, oNetwork.getEdges().size());
			}
			else
			{
				assertNoOverlapNetwork(output);
			}
		}
	}

	/**
	 * Private method that asserts an error will result from calling the oneMulti calculation
	 */
	private void assertException(Project p, CalculateDiff calc, int ID, boolean diff, boolean overlap)
	{
		boolean exceptionThrown = false;
		try
		{
			calc.calculateOneDifferentialNetwork(p, ID, cutoff, diff, overlap);
		}
		catch (IllegalArgumentException e)
		{
			exceptionThrown = true;
		}
		assertTrue(exceptionThrown);
	}

	/**
	 * Private method that asserts that a pair can not be queried from this output result
	 * @param output the differential output object
	 */
	private void assertNoPair(DifferentialOutput output)
	{
		boolean exceptionThrown = false;
		try
		{
			@SuppressWarnings("unused")
			OutputNetworkPair pair = output.getOutputAsPair();
		}
		catch (IllegalArgumentException e)
		{
			exceptionThrown = true;
		}
		assertTrue(exceptionThrown);
	}

	/**
	 * Private method that asserts that a differential network can not be queried from this output result
	 * @param output the differential output object
	 */
	private void assertNoDiffNetwork(DifferentialOutput output)
	{
		boolean exceptionThrown = false;
		try
		{
			@SuppressWarnings("unused")
			DifferentialNetwork dn = output.getDifferentialNetwork();
		}
		catch (IllegalArgumentException e)
		{
			exceptionThrown = true;
		}
		assertTrue(exceptionThrown);
	}

	/**
	 * Private method that asserts that an overlap network can not be queried from this output result
	 * @param output the differential output object
	 */
	private void assertNoOverlapNetwork(DifferentialOutput output)
	{
		boolean exceptionThrown = false;
		try
		{
			@SuppressWarnings("unused")
			OverlappingNetwork dn = output.getOverlappingNetwork();
		}
		catch (IllegalArgumentException e)
		{
			exceptionThrown = true;
		}
		assertTrue(exceptionThrown);
	}

}
