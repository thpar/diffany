package be.svlandeg.diffany.junit;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.examples.FuzzyNetworks;


/**
 * Class that automatically tests the outputs of some examples of 'fuzzy' differential networks.
 * 
 * @author Sofie Van Landeghem
 */
public class TestFuzzyDiff extends TestGeneric
{
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct differential results with fuzziness overlap factor 3 out of 3.
	 */
	@Test
	public void testFuzzyDiff_4()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 4;
		int ID = ex.getTestConfigurationWithReference(p, overlap_cutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, 70, -1, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 0);
		assertNrDiffNetworks(output, 1);

		// Testing the edges in the overlap network
		DifferentialNetwork dn = output.getDifferentialNetworks().iterator().next();

		Set<Edge> sEdges = dn.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(dn, "X", "Y", false, "increases_unspecified_regulation", false, 0.1);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct differential results with fuzziness overlap factor 2 out of 3.
	 */
	@Test
	public void testFuzzyDiff_3()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 3;
		int ID = ex.getTestConfigurationWithReference(p, overlap_cutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, 75, -1, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 0);
		assertNrDiffNetworks(output, 1);

		// Testing the edges in the overlap network
		DifferentialNetwork dn = output.getDifferentialNetworks().iterator().next();

		Set<Edge> sEdges = dn.getEdges();
		assertEquals(2, sEdges.size());

		assertAnEdge(dn, "B", "A", false, "increases_ppi", false, 0.4);
		assertAnEdge(dn, "X", "Y", false, "increases_regulation", false, 1.2);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct differential results with fuzziness overlap factor 2 out of 3.
	 */
	@Test
	public void testFuzzyDiff_2()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getProject();
		int overlap_cutoff = 2;
		int ID = ex.getTestConfigurationWithReference(p, overlap_cutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, 75, -1, true);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrOverlapNetworks(output, 0);
		assertNrDiffNetworks(output, 1);

		// Testing the edges in the differential network
		DifferentialNetwork dn = output.getDifferentialNetworks().iterator().next();

		Set<Edge> sEdges = dn.getEdges();
		assertEquals(4, sEdges.size());

		assertAnEdge(dn, "B", "A", false, "increases_ppi", false, 0.8);
		assertAnEdge(dn, "A", "B", false, "increases_colocalization", false, 0.2); 
		
		assertAnEdge(dn, "X", "Y", false, "increases_regulation", false, 1.3);
		
		assertAnEdge(dn, "M", "N", false, "increases_phosphorylation", false, 0.7);
	}
	
}
