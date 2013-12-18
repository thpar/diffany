package be.svlandeg.diffany.cytoscape.vizmapper;

import java.awt.Color;
import java.awt.Paint;
import java.util.HashSet;
import java.util.Set;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.view.model.DiscreteRange;
import org.cytoscape.view.model.Range;
import org.cytoscape.view.presentation.property.BasicVisualLexicon;
import org.cytoscape.view.presentation.property.NodeShapeVisualProperty;
import org.cytoscape.view.presentation.property.PaintVisualProperty;
import org.cytoscape.view.vizmap.VisualMappingFunctionFactory;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.view.vizmap.mappings.PassthroughMapping;

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
		this.specificMappings();
		
		services.getVisualMappingManager().addVisualStyle(this.vis);
	}

	private void defaultStyle() {		
		//network default style
		vis.setDefaultValue(BasicVisualLexicon.NETWORK_BACKGROUND_PAINT, Color.GRAY);
		
		//node default style
		vis.setDefaultValue(BasicVisualLexicon.NODE_SHAPE, NodeShapeVisualProperty.ROUND_RECTANGLE);
		vis.setDefaultValue(BasicVisualLexicon.NODE_FILL_COLOR, Color.YELLOW);
		
		//node basic mappings
		VisualMappingFunctionFactory vmffP = services.getVisualMappingFunctionFactory("passthrough");
		PassthroughMapping<String, ?> mapping = (PassthroughMapping<String, ?>)vmffP.createVisualMappingFunction(CyNetwork.NAME, String.class, BasicVisualLexicon.NODE_LABEL);
		vis.addVisualMappingFunction(mapping);	
		
		//edge style
		
	}

	/**
	 * Initialize this visual style
	 */
	abstract protected void specificMappings();
	
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
