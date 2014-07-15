package be.svlandeg.diffany.junit;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.OverlappingNetwork;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.examples.FuzzyOverlap;

/**
 * Class that automatically tests the outputs of some examples of 'fuzzy' overlap networks.
 * 
 * @author Sofie Van Landeghem
 */
public class TestFuzzyOverlap extends TestExamples
{

	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 4 out of 4 100%.
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyOverlapWithoutReference_4()
	{
		FuzzyOverlap ex = new FuzzyOverlap();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 4;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff);
		
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, false, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		for (Edge e : sEdges)
		{
			System.out.println(e);
		}
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "X", "Y", false, false, "regulation", false, 0.3);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 3 out of 4 (75%).
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyOverlapWithoutReference_3()
	{
		FuzzyOverlap ex = new FuzzyOverlap();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 3;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff);
		
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, false, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		for (Edge e : sEdges)
		{
			System.out.println(e);
		}
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "A", "B", false, false, "ppi", false, 0.3);
		assertAnEdge(on, "X", "Y", false, false, "regulation", false, 0.3);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 2 out of 4 (50%).
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */
	@Test
	public void testFuzzyOverlapWithoutReference_2()
	{
		FuzzyOverlap ex = new FuzzyOverlap();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 3;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff);
		
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, false, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		for (Edge e : sEdges)
		{
			System.out.println(e);
		}
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "A", "B", false, false, "ppi", false, 0.6);
		assertAnEdge(on, "B", "A", false, false, "ppi", false, 0.4);
		
		//assertAnEdge(on, "X", "Y", false, false, "positive regulation", false, 0.3);	will be removed after clean-up of the network because 'regulation' has higher weight
		assertAnEdge(on, "X", "Y", false, false, "regulation", false, 0.6);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results with fuzziness overlap factor 1 out of 4 (25%).
	 * This method defines all input networks to be generic, i.e. there is not one specific reference network.
	 */ 
	@Test
	public void testFuzzyOverlapWithoutReference_1()
	{
		FuzzyOverlap ex = new FuzzyOverlap();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 3;
		int ID = ex.getTestConfigurationWithoutReference(p, overlap_cutoff);
		
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, false, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		for (Edge e : sEdges)
		{
			System.out.println(e);
		}
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "A", "B", false, false, "ppi", false, 0.8);
		assertAnEdge(on, "B", "A", false, false, "ppi", false, 0.8);
		
		assertAnEdge(on, "X", "Y", false, false, "positive regulation", false, 0.8);
		assertAnEdge(on, "X", "Y", false, false, "negative regulation", false, 0.5);
	}
	
	/**
	 * TODO: different factors
	 * JUNIT Test: check whether the example network with multiple conditions produces correct results when varying the fuzziness overlap factor.
	 * This method defines one of the networks to be a reference network.
	 */
	@Test
	public void testFuzzyOverlapWithReference()
	{
		FuzzyOverlap ex = new FuzzyOverlap();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 4;
		int ID = ex.getTestConfigurationWithReference(p, overlap_cutoff);
		
		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, false, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 1);

		// Testing the edges in the overlap network
		OverlappingNetwork on = output.getOverlappingNetworks().iterator().next();

		Set<Edge> sEdges = on.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(on, "X", "Y", false, false, "regulation", false, 0.3);
	}
}
