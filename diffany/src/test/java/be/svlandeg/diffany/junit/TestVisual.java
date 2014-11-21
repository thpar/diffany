package be.svlandeg.diffany.junit;

import static org.junit.Assert.assertEquals;

import java.awt.Color;

import org.junit.Test;

import be.svlandeg.diffany.core.networks.Condition;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;
import be.svlandeg.diffany.core.visualstyle.EdgeStyle.ArrowHead;

/** 
 * Class that automatically tests the visualisation properties of the DefaultEdgeOntology.
 * 
 * @author Sofie Van Landeghem
 */
public class TestVisual
{

	/**
	 * Test the assigned colors and arrowheads 
	 * in the source and differential visualisation styles
	 */
	@Test
	public void testAll()
	{
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project("testProject", eo);
		
		// ****** SOURCE NETWORK ****** //
		
		ReferenceNetwork ref = new ReferenceNetwork("testRef", 1, null);
		ref.addEdge(new Edge("somethingRandom", new Node("A", "A"), new Node("B", "B"), true));
		ref.addEdge(new Edge("increase", new Node("C", "C"), new Node("D", "D"), true));
		ref.addEdge(new Edge("ppiptm", new Node("E", "E"), new Node("F", "F"), true));
		ref.addEdge(new Edge("increasewhatever", new Node("G", "G"), new Node("H", "H"), true));
		
		Condition c = new Condition("sample condition");
		ConditionNetwork cn = new ConditionNetwork("testCon", 2, null, c);
		cn.addEdge(new Edge("increase_ppi", new Node("I", "I"), new Node("J", "J"), true));
		cn.addEdge(new Edge("neutral", new Node("K", "K"), new Node("L", "L"), true));
		
		boolean cleanInput = true;
		p.addRunConfiguration(ref, cn, cleanInput, null);

		// process types
		assertColorInSource("ppi" , Color.YELLOW, eo);
		assertColorInSource("ptm" , Color.BLUE, eo);
		assertColorInSource("phosphorylates" , Color.BLUE, eo);
		assertColorInSource("ubiquitinate" , Color.BLUE, eo);
		assertColorInSource("methylation" , Color.BLUE, eo);
		
		assertArrowHeadInSource("ppi" , ArrowHead.NONE, eo);
		assertArrowHeadInSource("ptm" , ArrowHead.DIAMOND, eo);
		assertArrowHeadInSource("phosphorylates" , ArrowHead.DIAMOND, eo);
		assertArrowHeadInSource("ubiquitinate" , ArrowHead.DIAMOND, eo);
		assertArrowHeadInSource("methylation" , ArrowHead.DIAMOND, eo);
		
		// activity flow types
		assertColorInSource("negative gi" , Color.RED, eo);
		assertColorInSource("negative regulation" , Color.RED, eo);
		assertColorInSource("negatively regulates" , Color.RED, eo);
		assertColorInSource("negative" , Color.RED, eo);
		assertColorInSource("pos-reg", Color.GREEN, eo);
		assertColorInSource("positive regulation", Color.GREEN, eo);
		
		assertArrowHeadInSource("negative gi" , ArrowHead.NONE, eo);
		assertArrowHeadInSource("negative regulation" , ArrowHead.T, eo);
		assertArrowHeadInSource("negatively regulates" , ArrowHead.T, eo);
		assertArrowHeadInSource("negative" , ArrowHead.T, eo);
		assertArrowHeadInSource("pos-reg", ArrowHead.ARROW, eo);
		assertArrowHeadInSource("positive regulation", ArrowHead.ARROW, eo);
		
		// unmapped or neutral types
		assertColorInSource("regulates" , Color.LIGHT_GRAY, eo);
		assertColorInSource("regulate", Color.LIGHT_GRAY, eo);

		assertArrowHeadInSource("regulates" , ArrowHead.ARROW, eo);
		assertArrowHeadInSource("regulate", ArrowHead.ARROW, eo);
		
		assertColorInSource("somethingRandom", Color.LIGHT_GRAY, eo);
		assertColorInSource("increase", Color.LIGHT_GRAY, eo);
		assertColorInSource("ppiptm", Color.LIGHT_GRAY, eo);
		assertColorInSource("increasewhatever", Color.LIGHT_GRAY, eo);
		assertColorInSource("increase_ppi", Color.LIGHT_GRAY, eo);
		assertColorInSource("somethingrandom", Color.LIGHT_GRAY, eo);
		assertColorInSource("neutral", Color.LIGHT_GRAY, eo);
		
		assertArrowHeadInSource("somethingRandom", ArrowHead.NONE, eo);
		assertArrowHeadInSource("increase", ArrowHead.NONE, eo);
		assertArrowHeadInSource("ppiptm", ArrowHead.NONE, eo);
		assertArrowHeadInSource("increasewhatever", ArrowHead.NONE, eo);
		assertArrowHeadInSource("increase_ppi", ArrowHead.NONE, eo);
		assertArrowHeadInSource("somethingrandom", ArrowHead.NONE, eo);
		assertArrowHeadInSource("neutral", ArrowHead.NONE, eo);
		
		// ****** DIFFERENTIAL NETWORK ****** //
		
		// process types
		assertColorInDifferential("increase_ppi", new Color(51,255,0), eo); // cf. DefaultEdgeOntology.addDefaultPaintedParents
		assertColorInDifferential("increase_whatever", Color.GREEN, eo);
		assertColorInDifferential("decrease_ptm", Color.RED, eo);
		assertColorInDifferential("decreases_ptm", Color.RED, eo);
		
		assertArrowHeadInDifferential("increase_ppi", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("increase_whatever", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("decrease_ptm", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("decreases_ptm", ArrowHead.ARROW, eo);
		
		// activity flow types
		
		assertColorInDifferential("decreases_regulation",  new Color(204,0,0), eo);	// cf. DefaultEdgeOntology.addDefaultPaintedParents
		assertColorInDifferential("decrease_regulation", Color.RED, eo);
		assertColorInDifferential("increases_regulation", new Color(0,153,0), eo);	// cf. DefaultEdgeOntology.addDefaultPaintedParents
		assertColorInDifferential("increase_regulation", Color.GREEN, eo);
		
		assertArrowHeadInDifferential("decreases_regulation", ArrowHead.ARROW, eo);
		assertArrowHeadInDifferential("decrease_regulation", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("increases_regulation", ArrowHead.ARROW, eo);
		assertArrowHeadInDifferential("increase_regulation", ArrowHead.NONE, eo);
		
		// unmapped types
		assertColorInDifferential("increase", Color.GRAY, eo);
		assertColorInDifferential("decrease", Color.GRAY, eo);
		assertColorInDifferential("increases", Color.GRAY, eo);
		assertColorInDifferential("decreases", Color.GRAY, eo);
		assertArrowHeadInDifferential("increase", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("decrease", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("increases", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("decreases", ArrowHead.NONE, eo);
		
		assertColorInDifferential("regulates", Color.GRAY, eo);
		assertColorInDifferential("somethingRandom", Color.GRAY, eo);
		assertColorInDifferential("positive", Color.GRAY, eo);
		assertColorInDifferential("negativeSomething", Color.GRAY, eo);
		assertColorInDifferential("ppiptm", Color.GRAY, eo);
		assertColorInDifferential("somethingrandom", Color.GRAY, eo);
		assertColorInDifferential("ppi", Color.GRAY, eo);
		assertColorInDifferential("ptm", Color.GRAY, eo);
		assertColorInDifferential("regulate", Color.GRAY, eo);
		assertColorInDifferential("neutral", Color.GRAY, eo);
		
		assertArrowHeadInDifferential("regulates", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("somethingRandom", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("positive", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("negativeSomething", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("ppiptm", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("somethingrandom", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("ppi", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("ptm", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("regulate", ArrowHead.NONE, eo);
		assertArrowHeadInDifferential("neutral", ArrowHead.NONE, eo);
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
		Color p = eo.getSourceEdgeDrawing().getEdgeStyle(type).getColor();
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
		Color p = eo.getDifferentialEdgeDrawing().getEdgeStyle(cat).getColor();
		assertEquals(c, p);
	}
	
	/**
	 * Private method that asserts whether a certain interaction equals a certain arrowhead in a source network.
	 * 
	 * @param type the interaction type in the source network
	 * @param a the arrowhead the edge should get
	 * @param eo the edge ontology that will derive the color
	 */
	private void assertArrowHeadInSource(String type, ArrowHead a, EdgeOntology eo)
	{
		ArrowHead p = eo.getSourceEdgeDrawing().getEdgeStyle(type).getArrowHead();
		assertEquals(p, a);
	}
	
	/**
	 * Private method that asserts whether a certain interaction equals a certain arrowhead in a differential network.
	 * 
	 * @param cat the interaction cat in the differential network
	 * @param a the arrowhead the edge should get
	 * @param eo the edge ontology that will derive the color
	 */
	private void assertArrowHeadInDifferential(String cat, ArrowHead a, EdgeOntology eo)
	{
		ArrowHead p = eo.getDifferentialEdgeDrawing().getEdgeStyle(cat).getArrowHead();
		assertEquals(a, p);
	}
	
}
