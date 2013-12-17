package be.svlandeg.diffany.cytoscape.vizmapper;

import java.awt.Color;
import java.awt.Paint;
import java.util.HashSet;
import java.util.Set;

import org.cytoscape.model.CyNode;
import org.cytoscape.view.model.DiscreteRange;
import org.cytoscape.view.model.Range;
import org.cytoscape.view.presentation.property.NodeShapeVisualProperty;
import org.cytoscape.view.presentation.property.PaintVisualProperty;
import org.cytoscape.view.vizmap.VisualStyle;

import be.svlandeg.diffany.internal.Services;

public abstract class AbstractVisualDiffanyStyle {

	protected String name;
	protected VisualStyle vis;
	protected Services services;
	
	/**
	 * Create a new visual style and initialize it according to its type. Then register the
	 * new style with Cytoscape.
	 * 
	 * @param name the displayed name in the VizMapper
	 * @param services the services object
	 */
	public AbstractVisualDiffanyStyle(String name, Services services){
		this.name = name;
		this.services = services;
		this.vis = services.getVisualStyleFactory().createVisualStyle(name);
		
		this.defaultStyle();
		this.init();
		
		services.getVisualMappingManager().addVisualStyle(this.vis);
	}

	private void defaultStyle() {
		System.out.println("Added vis defaults");
		
		NodeShapeVisualProperty nodeShape = 
				new NodeShapeVisualProperty(NodeShapeVisualProperty.HEXAGON, "shape", "shaaaape", CyNode.class);
		vis.setDefaultValue(nodeShape, NodeShapeVisualProperty.HEXAGON);
		
		Set<Paint> colorSet = new HashSet<Paint>();
		colorSet.add(Color.GREEN);
		colorSet.add(Color.YELLOW);
		Range<Paint> range = new DiscreteRange<Paint>(Paint.class, colorSet);
		
		PaintVisualProperty nodePaint = 
				new PaintVisualProperty(Color.GREEN, range, "color", "colooooor", CyNode.class);
		vis.setDefaultValue(nodePaint, Color.GREEN);
		
		
	}

	/**
	 * Initialize this visual style
	 */
	abstract protected void init();
	
	/**
	 * Re-initialize mappings according to the content of selected networks. 
	 */
	abstract public void refresh(); 
	
	/**
	 * 
	 * @return the actual {@link VisualStyle} object.
	 */
	public VisualStyle getVisualStyle(){
		return this.vis;
	}
	
	private Set<String> getAllInteractions(){
		Set<String> interactions = new HashSet<String>();
		
		return interactions;
	}
}
