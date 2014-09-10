package be.svlandeg.diffany.junit;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.examples.OsmoticSampleTest;

/** 
 * This class provides examples to benchmark the fuzzy consensus & differential functionality, using real examples from the osmotic data.
 * 
 * @author Sofie Van Landeghem
 */
public class TestOsmoticSample extends TestGeneric
{
	 
	/**
	 * JUNIT Test: check whether the example fuzzy network produces correct results with the min operator and cutoff 4 out of 4
	 */
	@Test
	public void testFuzzy_4_min()
	{
		OsmoticSampleTest ex = new OsmoticSampleTest();
		double weight_cutoff = 0.0;
		Project p = ex.getTestProject();
		int supportingCutoff = 5;
		int ID = ex.getTestDiffConfiguration(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, 70, 80, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);
		
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();

		// Testing the edges in the differential network
		DifferentialNetwork dn = pair.getDifferentialNetwork();
		Set<Edge> dEdges = dn.getEdges();
		assertEquals(2, dEdges.size());
		
		assertAnEdge(dn, "M", "N", false, "decreases_regulation", false, 1.1);
		assertAnEdge(dn, "O", "P", false, "decreases_regulation", false, 1.1);
	}

}
