package be.svlandeg.diffany.junit;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import org.junit.Test;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.networks.DifferentialNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.OutputNetworkPair;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunOutput;
import be.svlandeg.diffany.examples.FuzzyNetworks2;

/** 
 * This class provides examples to benchmark the fuzzy consensus & differential functionality, varying consensus cutoffs and min/max operators.  
 * 
 * @author Sofie Van Landeghem
 */
public class TestFuzzyBoth extends TestGeneric
{

	/**
	 * JUNIT Test: check whether the example fuzzy network produces correct results with the min operator and cutoff 4 out of 4
	 */
	@Test
	public void testFuzzy_4_min()
	{
		FuzzyNetworks2 ex = new FuzzyNetworks2();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int supportingCutoff = 4;
		int ID = ex.getDefaultRunConfigurationID(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, 70, 80, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);
		
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();

		// Testing the edges in the differential network
		DifferentialNetwork dn = pair.getDifferentialNetwork();
		Set<Edge> dEdges = dn.getEdges();
		assertEquals(1, dEdges.size());
		assertAnEdge(dn, "X", "Y", false, "decreases_regulation", false, 5); 
		
		// Testing the edges in the consensus network
		ConsensusNetwork cn = pair.getConsensusNetwork();
		Set<Edge> cEdges = cn.getEdges();
		assertEquals(1, cEdges.size());
		assertAnEdge(cn, "A", "B", true, "ppi", false, 0.4); 
	}
	
	/**
	 * JUNIT Test: check whether the example fuzzy network produces correct results with the min operator and cutoff 4 out of 4
	 */
	@Test
	public void testFuzzy_4_max()
	{
		FuzzyNetworks2 ex = new FuzzyNetworks2();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int supportingCutoff = 4;
		int ID = ex.getDefaultRunConfigurationID(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, 70, 80, false, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);
		
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();

		// Testing the edges in the differential network
		DifferentialNetwork dn = pair.getDifferentialNetwork();
		Set<Edge> dEdges = dn.getEdges();
		assertEquals(1, dEdges.size());
		assertAnEdge(dn, "X", "Y", false, "decreases_regulation", false, 5); 
		
		// Testing the edges in the consensus network
		ConsensusNetwork cn = pair.getConsensusNetwork();
		Set<Edge> cEdges = cn.getEdges();
		assertEquals(1, cEdges.size());
		assertAnEdge(cn, "A", "B", true, "ppi", false, 1.2); 
	}
	
	/**
	 * JUNIT Test: check whether the example fuzzy network produces correct results with the min operator and cutoff 3 out of 4
	 */
	@Test
	public void testFuzzy_3_min()
	{
		FuzzyNetworks2 ex = new FuzzyNetworks2();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int supportingCutoff = 3;
		int ID = ex.getDefaultRunConfigurationID(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, 70, 80, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);
		
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();

		// Testing the edges in the differential network
		DifferentialNetwork dn = pair.getDifferentialNetwork();
		Set<Edge> dEdges = dn.getEdges();
		assertEquals(2, dEdges.size());
		assertAnEdge(dn, "A", "B", true, "decrease_ppi", false, 0.2); 
		assertAnEdge(dn, "X", "Y", false, "decreases_regulation", false, 5); 
		
		// Testing the edges in the consensus network
		ConsensusNetwork cn = pair.getConsensusNetwork();
		Set<Edge> cEdges = cn.getEdges();
		assertEquals(2, cEdges.size());
		assertAnEdge(cn, "A", "B", true, "ppi", false, 0.6); 
		assertAnEdge(cn, "X", "Y", false, "regulation", false, 4); 
	}
	
	/**
	 * JUNIT Test: check whether the example fuzzy network produces correct results with the min operator and cutoff 3 out of 4
	 */
	@Test
	public void testFuzzy_3_max()
	{
		FuzzyNetworks2 ex = new FuzzyNetworks2();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int supportingCutoff = 3;
		int ID = ex.getDefaultRunConfigurationID(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, 70, 80, false, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);
		
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();

		// Testing the edges in the differential network
		DifferentialNetwork dn = pair.getDifferentialNetwork();
		Set<Edge> dEdges = dn.getEdges();
		assertEquals(2, dEdges.size());
		assertAnEdge(dn, "A", "B", true, "decrease_ppi", false, 0.2);
		assertAnEdge(dn, "X", "Y", false, "decreases_regulation", false, 5); 
		
		// Testing the edges in the consensus network
		ConsensusNetwork cn = pair.getConsensusNetwork();
		Set<Edge> cEdges = cn.getEdges();
		assertEquals(2, cEdges.size());
		assertAnEdge(cn, "A", "B", true, "ppi", false, 1.2); 
		assertAnEdge(cn, "X", "Y", false, "regulation", false, 11); 
	}
	
	/**
	 * JUNIT Test: check whether the example fuzzy network produces correct results with the min operator and cutoff 2 out of 4
	 */
	@Test
	public void testFuzzy_2_min()
	{
		FuzzyNetworks2 ex = new FuzzyNetworks2();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int supportingCutoff = 2;
		int ID = ex.getDefaultRunConfigurationID(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, 70, 80, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);
		
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();

		// Testing the edges in the differential network
		DifferentialNetwork dn = pair.getDifferentialNetwork();
		Set<Edge> dEdges = dn.getEdges();
		assertEquals(1, dEdges.size());
		assertAnEdge(dn, "X", "Y", false, "decreases_regulation", false, 5); 
		
		// Testing the edges in the consensus network
		ConsensusNetwork cn = pair.getConsensusNetwork();
		Set<Edge> cEdges = cn.getEdges();
		assertEquals(2, cEdges.size());
		assertAnEdge(cn, "A", "B", true, "ppi", false, 0.8); 
		assertAnEdge(cn, "X", "Y", false, "regulation", false, 6); 
	}
	
	/**
	 * JUNIT Test: check whether the example fuzzy network produces correct results with the min operator and cutoff 3 out of 4
	 */
	@Test
	public void testFuzzy_2_max()
	{
		FuzzyNetworks2 ex = new FuzzyNetworks2();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int supportingCutoff = 2;
		int ID = ex.getDefaultRunConfigurationID(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, 70, 80, false, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrPairs(output, 1);
		
		OutputNetworkPair pair = output.getOutputAsPairs().iterator().next();

		// Testing the edges in the differential network
		DifferentialNetwork dn = pair.getDifferentialNetwork();
		Set<Edge> dEdges = dn.getEdges();
		assertEquals(1, dEdges.size());
		assertAnEdge(dn, "X", "Y", false, "decreases_regulation", false, 5); 
		
		// Testing the edges in the consensus network
		ConsensusNetwork cn = pair.getConsensusNetwork();
		Set<Edge> oEdges = cn.getEdges();
		assertEquals(2, oEdges.size());
		assertAnEdge(cn, "A", "B", true, "ppi", false, 1.2); 
		assertAnEdge(cn, "X", "Y", false, "regulation", false, 11); 
	}
}
