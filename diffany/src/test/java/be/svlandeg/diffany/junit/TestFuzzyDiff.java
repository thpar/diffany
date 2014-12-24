package be.svlandeg.diffany.junit;

/*
 * #%L
 * Diffany
 * %%
 * Copyright (C) 2014 PSB/UGent - Sofie Van Landeghem and Thomas Van Parys
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */


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
	 * JUNIT Test: check whether the example network with multiple conditions produces correct differential results with fuzziness supporting networks factor 3 out of 3.
	 */
	@Test
	public void testFuzzyDiff_4()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int supportingCutoff = 4;
		int ID = ex.getTestConfigurationWithReference(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, null, null, 70, -1, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 0);
		assertNrDiffNetworks(output, 1);

		// Testing the edges in the differential network
		DifferentialNetwork dn = output.getDifferentialNetworks().iterator().next();

		Set<Edge> sEdges = dn.getEdges();
		assertEquals(1, sEdges.size());

		assertAnEdge(dn, "X", "Y", false, "increases_unspecified_regulation", false, 0.1);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct differential results with fuzziness supporting networks factor 2 out of 3.
	 */
	@Test
	public void testFuzzyDiff_3()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int supportingCutoff = 3;
		int ID = ex.getTestConfigurationWithReference(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, null, null, 75, -1, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 0);
		assertNrDiffNetworks(output, 1);

		// Testing the edges in the differential network
		DifferentialNetwork dn = output.getDifferentialNetworks().iterator().next();

		Set<Edge> sEdges = dn.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(dn, "A", "B", false, "decreases_colocalization", false, 0.3);
		assertAnEdge(dn, "B", "A", false, "increases_ppi", false, 0.4);
		assertAnEdge(dn, "X", "Y", false, "increases_regulation", false, 1.2);
	}
	
	/**
	 * JUNIT Test: check whether the example network with multiple conditions produces correct differential results with fuzziness supporting networks factor 2 out of 3.
	 */
	@Test
	public void testFuzzyDiff_2()
	{
		FuzzyNetworks ex = new FuzzyNetworks();
		double weight_cutoff = 0.0;
		Project p = ex.getDefaultProject();
		int supportingCutoff = 2;
		int ID = ex.getTestConfigurationWithReference(p, supportingCutoff);

		new CalculateDiff().calculateOneDifferentialNetwork(p, ID, weight_cutoff, null, null, 75, -1, true, null);

		// Testing that there is exactly one differential network created
		RunOutput output = p.getOutput(ID);
		assertNrConsensusNetworks(output, 0);
		assertNrDiffNetworks(output, 1);

		// Testing the edges in the differential network
		DifferentialNetwork dn = output.getDifferentialNetworks().iterator().next();

		Set<Edge> sEdges = dn.getEdges();
		assertEquals(3, sEdges.size());

		assertAnEdge(dn, "B", "A", false, "increases_ppi", false, 0.8);
		assertAnEdge(dn, "X", "Y", false, "increases_regulation", false, 1.3);
		assertAnEdge(dn, "M", "N", false, "increases_phosphorylation", false, 0.7);
	}
	
}
