package be.svlandeg.diffany.junit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.Map;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.examples.MultipleConditionTest;

/**
 * Class that automatically tests the RunConfiguration settings and results stored in a Project entity.
 * This class only checks the type of networks generated, as other tests examine the content of those networks specifically.
 * 
 * @author Sofie Van Landeghem
 */
public class TestConfiguration extends TestGeneric
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

		// 1-all differential configuration : differential+consensus calculation
		testDifferentialOneMulti(p, calc, ex, calls, 10, 11);
		calls++;

		// (the exact same) 1-all differential+consensus calculation
		testDifferentialOneMulti(p, calc, ex, calls, 12, 13);
		calls++;

		// 1-all differential (only) calculation
		testDifferentialOneMulti(p, calc, ex, calls, 14, -1);
		calls++;

		// 1-all consensus (only) calculation - with differential configuration
		testDifferentialOneMulti(p, calc, ex, calls, -1, 15);
		calls++;

		// 1-all calculation which doesn't actually do anything (no diff, no consensus)
		testDifferentialOneMulti(p, calc, ex, calls, -1, -1);
		calls++;

		// 1-all consensus configuration - assess correct failure of differential calculation
		testConsensusOneMulti(p, calc, ex, calls, 20);
		calls++;

		// pairwise differential configuration : differential+consensus calculation
		testDifferentialPairwise(p, calc, ex, calls, true, true, 30);
		calls++;

		// pairwise differential (only) calculation
		testDifferentialPairwise(p, calc, ex, calls, true, false, 40);
		calls++;

		// pairwise consensus (only) calculation - with differential configuration
		testDifferentialPairwise(p, calc, ex, calls, false, false, 50);
		calls++;

		// pairwise calculation which doesn't actually do anything (no diff, no consensus)
		testDifferentialPairwise(p, calc, ex, calls, false, false, 60);
		calls++;
	}

	/**
	 * Test an initial empty project: should contain no differential output.
	 */
	private void testIni(Project p, MultipleConditionTest ex, int calls)
	{
		assertTrue(0 == p.getNumberOfRuns());

		int ID = ex.getTestDiffConfiguration(p);
		assertTrue(calls + 1 == p.getNumberOfRuns());
		
		RunOutput dOutput = p.getOutput(ID);
		assertNrDiffNetworks(dOutput, 0);
	}

	/**
	 * Test a 1-multi differential configuration, with differential and/or consensus calculation.
	 */
	private void testDifferentialOneMulti(Project p, CalculateDiff calc, MultipleConditionTest ex, int calls, int diffID, int consensusID)
	{
		int ID = ex.getTestDiffConfiguration(p);
		assertTrue(calls + 1 == p.getNumberOfRuns());

		// it should not make a difference how many times this method is called!
		calc.calculateOneDifferentialNetwork(p, ID, cutoff, diffID, consensusID, true, null);
		calc.calculateOneDifferentialNetwork(p, ID, cutoff, diffID, consensusID, true, null);
		RunOutput output = p.getOutput(ID);
		
		boolean diff = diffID >= 0;
		boolean consensus = consensusID >= 0;
		
		if (diff && consensus)
		{
			assertNrPairs(output, 1);
			OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();
			
			DifferentialNetwork dNetwork = pair.getDifferentialNetwork();
			assertEquals(8, dNetwork.getEdges().size());

			ConsensusNetwork oNetwork = pair.getConsensusNetwork();
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
		if (consensus)
		{
			assertNrConsensusNetworks(output, 1);
			assertEquals(1, output.getConsensusNetworks().size());
			ConsensusNetwork oNetwork = output.getConsensusNetworks().iterator().next();
			assertEquals(3, oNetwork.getEdges().size());
		}
		else
		{
			assertNrConsensusNetworks(output, 0);
		}
	}

	/**
	 * Test a 1-multi consensus configuration, should not contain differential networks.
	 */
	private void testConsensusOneMulti(Project p, CalculateDiff calc, MultipleConditionTest ex, int calls, int consensusID)
	{
		int ID = ex.getTestConsensusConfiguration(p);
		assertTrue(calls + 1 == p.getNumberOfRuns());

		assertException(p, calc, ID, consensusID++, consensusID++);
		assertException(p, calc, ID, consensusID++, -1);

		RunOutput dOutput = p.getOutput(ID);
		assertNrPairs(dOutput, 0);
		assertNrDiffNetworks(dOutput, 0);
		assertNrConsensusNetworks(dOutput, 0);

		calc.calculateOneDifferentialNetwork(p, ID, cutoff, -1, consensusID, true, null);
		dOutput = p.getOutput(ID);

		assertNrPairs(dOutput, 0);
		assertNrDiffNetworks(dOutput, 0);
		assertNrConsensusNetworks(dOutput, 1);

		ConsensusNetwork oNetwork = dOutput.getConsensusNetworks().iterator().next();
		assertEquals(3, oNetwork.getEdges().size());
	}

	/**
	 * Test a pairwise differential configuration, with differential and/or consensus calculation.
	 */
	private void testDifferentialPairwise(Project p, CalculateDiff calc, MultipleConditionTest ex, int calls, boolean diff, boolean consensus, int firstID)
	{
		int ID = ex.getTestDiffConfiguration(p);
		assertTrue(calls + 1 == p.getNumberOfRuns());

		calc.calculateAllPairwiseDifferentialNetworks(p, ID, cutoff, diff, consensus, firstID, true, null);
		RunOutput output = p.getOutput(ID);

		// First, test the number of result networks, depending on the settings.
		if (diff)
		{
			// ref vs. 2 conditions
			assertNrDiffNetworks(output, 2);
			if (consensus)
			{
				assertNrConsensusNetworks(output, 2);
			}
			else
			{
				assertNrConsensusNetworks(output, 0);
			}
		}
		
		else if (consensus)
		{
			// 3 pairwise comparisons
			assertNrPairs(output, 0);
			assertNrDiffNetworks(output, 0);
			assertNrConsensusNetworks(output, 3);
		}
		else
		{
			// no comparisons
			assertNrPairs(output, 0);
			assertNrDiffNetworks(output, 0);
			assertNrConsensusNetworks(output, 0);
		}

		// Now, test the number of edges per result network, depending on the settings.
		Map<String, Integer> diffcounts = new HashMap<String, Integer>();
		diffcounts.put("diff_Draughty", 12);
		diffcounts.put("diff_Salty", 11);

		Map<String, Integer> consensuscounts = new HashMap<String, Integer>();
		consensuscounts.put("consensus_Draughty_Reference", 3);
		consensuscounts.put("consensus_Reference_Salty", 5);
		consensuscounts.put("consensus_Draughty_Salty", 5);	

		if (diff && consensus)
		{
			for (OutputNetworkPair pair : output.getOutputAsPairs())
			{
				DifferentialNetwork dNetwork = pair.getDifferentialNetwork();
				int countD = diffcounts.get(dNetwork.getName());
				assertEquals(countD, dNetwork.getEdges().size());

				ConsensusNetwork oNetwork = pair.getConsensusNetwork();
				int countO = consensuscounts.get(oNetwork.getName());
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
		if (consensus)
		{
			for (ConsensusNetwork oNetwork : output.getConsensusNetworks())
			{
				int countO = consensuscounts.get(oNetwork.getName());
				assertEquals(countO, oNetwork.getEdges().size());
			}
		}
	}

	/**
	 * Private method that asserts an error will result from calling the oneMulti calculation with a specific runconfiguration ID.
	 */
	private void assertException(Project p, CalculateDiff calc, int ID, int diffID, int consensusID)
	{
		boolean exceptionThrown = false;
		try
		{
			calc.calculateOneDifferentialNetwork(p, ID, cutoff, diffID, consensusID, true, null);
		}
		catch (IllegalArgumentException e)
		{
			exceptionThrown = true;
		}
		assertTrue(exceptionThrown);
	}
}
