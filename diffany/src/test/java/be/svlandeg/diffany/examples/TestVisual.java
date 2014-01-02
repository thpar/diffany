package be.svlandeg.diffany.examples;

import static org.junit.Assert.assertEquals;

import java.awt.Color;

import org.junit.Test;

import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.EdgeOntology;

/** 
 * Class that automatically tests the visualisation properties of the DefaultEdgeOntology.
 * 
 * @author Sofie Van Landeghem
 */
public class TestVisual
{

	/**
	 * Test the assigned colors in the source and differential visualisation styles
	 */
	@Test
	public void testColors()
	{
		EdgeOntology eo = new DefaultEdgeOntology();
		
		// ****** SOURCE NETWORK ****** //
		
		// process types
		assertColorInSource("ppi" , Color.CYAN, eo);
		assertColorInSource("ptm" , Color.BLUE, eo);
		assertColorInSource("phosphorylates" , Color.BLUE, eo);
		assertColorInSource("ubiquitinate" , Color.BLUE, eo);
		assertColorInSource("methylation" , Color.BLUE, eo);
		
		// activity flow types
		assertColorInSource("negatively regulates" , Color.RED, eo);
		assertColorInSource("negative" , Color.RED, eo);
		assertColorInSource("pos-reg", Color.GREEN, eo);
		assertColorInSource("positive regulation", Color.GREEN, eo);
		
		// unmapped types
		assertColorInSource("regulates" , Color.LIGHT_GRAY, eo);
		assertColorInSource("somethingRandom", Color.LIGHT_GRAY, eo);
		assertColorInSource("increase", Color.LIGHT_GRAY, eo);
		assertColorInSource("ppiptm", Color.LIGHT_GRAY, eo);
		assertColorInSource("increasewhatever", Color.LIGHT_GRAY, eo);
		assertColorInSource("increase_ppi", Color.LIGHT_GRAY, eo);
		assertColorInSource("somethingrandom", Color.LIGHT_GRAY, eo);
		assertColorInSource("regulate", Color.LIGHT_GRAY, eo);
		assertColorInSource("neutral", Color.LIGHT_GRAY, eo);
		assertColorInSource(null, Color.LIGHT_GRAY, eo);
		
		// ****** DIFFERENTIAL NETWORK ****** //
		
		// process types
		assertColorInDifferential("increase_ppi", Color.YELLOW, eo);
		assertColorInDifferential("increase_whatever", Color.YELLOW, eo);
		assertColorInDifferential("decrease_ptm", Color.ORANGE, eo);
		
		// activity flow types
		assertColorInDifferential("increase", Color.GREEN, eo);
		assertColorInDifferential("decrease", Color.RED, eo);
		
		// unmapped types
		assertColorInDifferential("regulates", Color.GRAY, eo);
		assertColorInDifferential("somethingRandom", Color.GRAY, eo);
		assertColorInDifferential("positive", Color.GRAY, eo);
		assertColorInDifferential("negativeSomething", Color.GRAY, eo);
		assertColorInDifferential("ppiptm", Color.GRAY, eo);
		assertColorInDifferential("increasewhatever", Color.GRAY, eo);
		assertColorInDifferential("somethingrandom", Color.GRAY, eo);
		assertColorInDifferential("ppi", Color.GRAY, eo);
		assertColorInDifferential("ptm", Color.GRAY, eo);
		assertColorInDifferential("regulate", Color.GRAY, eo);
		assertColorInDifferential("neutral", Color.GRAY, eo);
		assertColorInDifferential(null, Color.GRAY, eo);
	}
	
	
	/**
	 * Private method that asserts whether a certain interaction equals a certain color in a source network.
	 * 
	 * @param type the interaction type in the source network
	 * @param c the color the edge should get
	 * @param eo the edge ontology that will derive the color
	 */
	private void assertColorInSource(String type, Color c, EdgeOntology eo)
	{
		Color p = eo.getSourceEdgeStyle(type).getColor();
		assertEquals(c, p);
	}
	
	/**
	 * Private method that asserts whether a certain interaction equals a certain color in a differential network.
	 * 
	 * @param cat the interaction cat in the differential network
	 * @param c the color the edge should get
	 * @param eo the edge ontology that will derive the color
	 */
	private void assertColorInDifferential(String cat, Color c, EdgeOntology eo)
	{
		Color p = eo.getDifferentialEdgeStyle(cat).getColor();
		assertEquals(c, p);
	}
	
}
