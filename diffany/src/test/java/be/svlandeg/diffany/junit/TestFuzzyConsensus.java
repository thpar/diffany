package be.svlandeg.diffany.junit;

import static org.junit.Assert.*;

import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.examples.FuzzyNetworks;

/**
 * Class that automatically tests the outputs of some examples of 'fuzzy' consensus networks.
 * 
 * @author Sofie Van Landeghem
 */
public class TestFuzzyConsensus extends TestGeneric
{
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness supporting networks factor 4 out of 4 100%.
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyConsensusWithReferenceAsCondition_4()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int supportingCutoff = 4;
		int ID = ex.getTestConfigurationWithoutReference(p, supportingCutoff, true);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 15, true, null);

		// Testing that there is exactly one consensus network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 1);
		assertNrDiffNetworks(output, 0);
		assertNrDiffNetworks(output, 0);

		// Testing the edges in the consensus network
		ConsensusNetwork on = output.getConsensusNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.5);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness supporting networks factor 3 out of 4 (75%).
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyConsensusWithReferenceAsCondition_3()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int supportingCutoff = 3;
		int ID = ex.getTestConfigurationWithoutReference(p, supportingCutoff, true);

		
		boolean exception = false;
		try
		{
			new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 10, true, null);		// this is an ID of an input network and should thus not work!
		}
		catch (IllegalArgumentException e)
		{
			exception = true;
		}
		assertTrue(exception);
		
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 20, true, null);	

		// Testing that there is exactly one consensus network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 1);
		assertNrDiffNetworks(output, 0);
		assertNrDiffNetworks(output, 0);

		// Testing the edges in the consensus network
		ConsensusNetwork on = output.getConsensusNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(on, "A", "B", false, "colocalization", false, 0.3);
		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.6);
		assertAnEdge(on, "M", "N", false, "phosphorylation", true, 0.3);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness supporting networks factor 3 out of 4 (75%).
	 * Instead of the minimum operator, the maximum is used.
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyConsensusWithReferenceAsCondition_3_maxOperator()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int supportingCutoff = 3;
		int ID = ex.getTestConfigurationWithoutReference(p, supportingCutoff, true);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 30, false, null);

		// Testing that there is exactly one consensus network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 1);
		assertNrDiffNetworks(output, 0);
		assertNrDiffNetworks(output, 0);

		// Testing the edges in the consensus network
		ConsensusNetwork on = output.getConsensusNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(on, "A", "B", false, "colocalization", false, 0.8);
		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.8);
		assertAnEdge(on, "M", "N", false, "phosphorylation", true, 0.6);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness supporting networks factor 2 out of 4 (50%).
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyConsensusWithReferenceAsCondition_2()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int supportingCutoff = 2;
		int ID = ex.getTestConfigurationWithoutReference(p, supportingCutoff, true);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 40, true, null);
		
		// Testing that there is exactly one consensus network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 1);
		assertNrDiffNetworks(output, 0);
		assertNrDiffNetworks(output, 0);

		// Testing the edges in the consensus network
		ConsensusNetwork on = output.getConsensusNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(4, sEdges.size());

		assertAnEdge(on, "A", "B", false, "colocalization", false, 0.6);
		assertAnEdge(on, "B", "A", false, "ppi", false, 0.4);

		assertAnEdge(on, "X", "Y", false, "positive_regulation", false, 0.7);	
		
		assertAnEdge(on, "M", "N", false, "ptm", true, 0.5);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness supporting networks factor 1 out of 4 (25%).
	 * This should throw an error, as te support factor should always be above 1 (i.e. at least 2).
	 */
	@Test
	public void testFuzzyConsensusWithReferenceAsCondition_1()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		Project p = ex.getProject();
		int supportingCutoff = 1;
		
		boolean exception = false;
		try 
		{
			ex.getTestConfigurationWithoutReference(p, supportingCutoff, true);
		}
		catch(IllegalArgumentException e)
		{
			exception = true;
		}
		assertTrue(exception);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness supporting networks factor 0 out of 4 (0%). 
	 * This should throw an error, as the supporting networks factor should always be above 1 (i.e. at least 2).
	 */
	@Test
	public void testFuzzyConsensusWithReferenceAsCondition_0()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		Project p = ex.getProject();
		int supportingCutoff = 0;

		boolean exception = false;
		try
		{
			ex.getTestConfigurationWithoutReference(p, supportingCutoff, true);
		}
		catch (IllegalArgumentException e)
		{
			exception = true;
		}
		assertTrue(exception);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness supporting networks factor 4.
	 * This should throw an error, as there are only 3 input networks.
	 */
	@Test
	public void testFuzzyConsensusWithoutReference_4()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int supportingCutoff = 4;
		int ID = ex.getTestConfigurationWithoutReference(p, supportingCutoff, false);

		boolean exception = false;
		try
		{
			new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 15, true, null);
		}
		catch (IllegalArgumentException e)
		{
			exception = true;
		}
		assertTrue(exception);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness supporting networks factor 3 out of 3 (100%).
	 * This method uses only the condition-specific networks, all cast as generic input networks.
	 */
	@Test
	public void testFuzzyConsensusWithoutReference_3()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int supportingCutoff = 3;
		int ID = ex.getTestConfigurationWithoutReference(p, supportingCutoff, false);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 20, true, null);	

		// Testing that there is exactly one consensus network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 1);
		assertNrDiffNetworks(output, 0);
		assertNrDiffNetworks(output, 0);

		// Testing the edges in the consensus network
		ConsensusNetwork on = output.getConsensusNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.6);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness supporting networks factor 3 out of 3 (100%).
	 * Instead of the minimum operator, the maximum is used.
	 * This method uses only the condition-specific networks, all cast as generic input networks.
	 */
	@Test
	public void testFuzzyConsensusWithoutReference_3_maxOperator()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int supportingCutoff = 3;
		int ID = ex.getTestConfigurationWithoutReference(p, supportingCutoff, false);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 30, false, null);

		// Testing that there is exactly one consensus network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 1);
		assertNrDiffNetworks(output, 0);
		assertNrDiffNetworks(output, 0);

		// Testing the edges in the consensus network
		ConsensusNetwork on = output.getConsensusNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.8);
	}

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness supporting networks factor 2 out of 3 (67%).
	 * This method uses only the condition-specific networks, all cast as generic input networks.
	 */
	@Test
	public void testFuzzyConsensusWithoutReference_2()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int supportingCutoff = 2;
		int ID = ex.getTestConfigurationWithoutReference(p, supportingCutoff, false);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 40, true, null);
		
		// Testing that there is exactly one consensus network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 1);
		assertNrDiffNetworks(output, 0);
		assertNrDiffNetworks(output, 0);

		// Testing the edges in the consensus network
		ConsensusNetwork on = output.getConsensusNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(4, sEdges.size());

		assertAnEdge(on, "A", "B", false, "ppi", false, 0.3);
		assertAnEdge(on, "B", "A", false, "ppi", false, 0.4);

		assertAnEdge(on, "X", "Y", false, "positive_regulation", false, 0.7);	
		
		assertAnEdge(on, "M", "N", false, "phosphorylation", true, 0.3);
	}


	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results when varying the fuzziness supporting networks factor.
	 * This method defines one of the networks to be a reference network.
	 */
	@Test
	public void testFuzzyConsensusWithReference_4()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int supportingCutoff = 4;
		int ID = ex.getTestConfigurationWithReference(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 70, true, null);

		// Testing that there is exactly one consensus network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 1);
		assertNrDiffNetworks(output, 0);
		assertNrDiffNetworks(output, 0);

		// Testing the edges in the consensus network
		ConsensusNetwork on = output.getConsensusNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.5);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results when varying the fuzziness supporting networks factor.
	 * This method defines one of the networks to be a reference network.
	 */
	@Test
	public void testFuzzyConsensusWithReference_3()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int supportingCutoff = 3;
		int ID = ex.getTestConfigurationWithReference(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 70, true, null);

		// Testing that there is exactly one consensus network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 1);
		assertNrDiffNetworks(output, 0);
		assertNrDiffNetworks(output, 0);

		// Testing the edges in the consensus network
		ConsensusNetwork on = output.getConsensusNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(on, "A", "B", false, "colocalization", false, 0.3);
		assertAnEdge(on, "X", "Y", false, "regulation", false, 0.5);		// important different with scenario without reference (0.6) !!!
		assertAnEdge(on, "M", "N", false, "phosphorylation", true, 0.3);
	}

	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results when varying the fuzziness supporting networks factor.
	 * This method defines one of the networks to be a reference network.
	 */
	@Test
	public void testFuzzyConsensusWithReference_2()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int supportingCutoff = 2;
		int ID = ex.getTestConfigurationWithReference(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, -1, 70, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 1);
		assertNrDiffNetworks(output, 0);
		assertNrDiffNetworks(output, 0);

		// Testing the edges in the consensus network
		ConsensusNetwork on = output.getConsensusNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(on, "A", "B", false, "colocalization", false, 0.6);
		assertAnEdge(on, "X", "Y", false, "negative_regulation", false, 0.5);		// important different with scenario without reference (pos 0.7) !!!
		assertAnEdge(on, "M", "N", false, "ptm", true, 0.5);
	}

	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results when varying the fuzziness supporting networks factor.
	 * This should throw an error, as the consensus factor should always be above 1 (i.e. at least 2).
	 */
	@Test
	public void testFuzzyConsensusWithReference_1()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		Project p = ex.getProject();
		int supportingCutoff = 1;
		boolean exception = false;
		try
		{
			ex.getTestConfigurationWithReference(p, supportingCutoff);
		}
		catch (IllegalArgumentException e)
		{
			exception = true;
		}
		assertTrue(exception);
	}
}
