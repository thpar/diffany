package be.svlandeg.diffany.cytoscape.vizmapper;

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

import java.awt.Color;
import java.awt.Paint;
import java.util.Set;

import org.cytoscape.model.CyEdge;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.presentation.property.BasicVisualLexicon;
import org.cytoscape.view.presentation.property.LineTypeVisualProperty;
import org.cytoscape.view.presentation.property.NodeShapeVisualProperty;
import org.cytoscape.view.presentation.property.values.ArrowShape;
import org.cytoscape.view.presentation.property.values.LineType;
import org.cytoscape.view.vizmap.VisualMappingFunction;
import org.cytoscape.view.vizmap.VisualMappingFunctionFactory;
import org.cytoscape.view.vizmap.VisualPropertyDependency;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.view.vizmap.mappings.DiscreteMapping;
import org.cytoscape.view.vizmap.mappings.PassthroughMapping;

import be.svlandeg.diffany.core.semantics.EdgeOntology;
import be.svlandeg.diffany.cytoscape.CyNetworkBridge;
import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.internal.Services;

/**
 * A {@link VisualStyle} wrapper with default values for all Diffany styles. 
 * All new styles within this app should extend this class.
 * 
 * @author Thomas Van Parys
 *
 */
public abstract class AbstractVisualDiffanyStyle {

	/**
	 * Name of the visual style
	 */
	protected String name;
	
	/**
	 * The {@link VisualStyle} itself.
	 */
	protected VisualStyle vis;
	
	/**
	 * App services
	 */
	protected Services services;
	
	//color definitions
	private static final Color NETWORK_BACKGROUND_COLOR = Color.WHITE;
	private static final Color NODE_COLOR = Color.GRAY;
	
	private static final double NODE_HEIGHT = 55;
	private static final double NODE_WIDTH = 80;
	private static final int FONT_SIZE = 11;
	
	
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

		registerStyle();
	}

	/**
	 * Registers this style with the VizMapper
	 */
	public void registerStyle(){
		services.getVisualMappingManager().addVisualStyle(this.vis);		
	}
	
	/**
	 * The default Diffany style that applies to all networks
	 */
	private void defaultStyle() {		
		//network default style
		vis.setDefaultValue(BasicVisualLexicon.NETWORK_BACKGROUND_PAINT, NETWORK_BACKGROUND_COLOR);
		
		//node default style
		vis.setDefaultValue(BasicVisualLexicon.NODE_SHAPE, NodeShapeVisualProperty.ELLIPSE);
		vis.setDefaultValue(BasicVisualLexicon.NODE_FILL_COLOR, NODE_COLOR);
		vis.setDefaultValue(BasicVisualLexicon.NODE_SELECTED_PAINT, NODE_COLOR.darker().darker());
		vis.setDefaultValue(BasicVisualLexicon.NODE_SIZE, 30d);
		vis.setDefaultValue(BasicVisualLexicon.NODE_BORDER_WIDTH, 0d);
		
		//Disable "Lock node width and height", so we can set a custom width and height.
		for(VisualPropertyDependency<?> visualPropertyDependency : vis.getAllVisualPropertyDependencies()) {
			if(visualPropertyDependency.getIdString().equals("nodeSizeLocked")) {
				visualPropertyDependency.setDependency(false);
				break;
			}
		}
		
		vis.setDefaultValue(BasicVisualLexicon.NODE_HEIGHT, NODE_HEIGHT);
		vis.setDefaultValue(BasicVisualLexicon.NODE_WIDTH, NODE_WIDTH);
		vis.setDefaultValue(BasicVisualLexicon.NODE_LABEL_FONT_SIZE, FONT_SIZE);
		
		
		
		//node basic mappings
		VisualMappingFunctionFactory vmffP = services.getVisualMappingFunctionFactory("passthrough");
		PassthroughMapping<String, ?> nodeLabelMapping = (PassthroughMapping<String, ?>)vmffP.createVisualMappingFunction(CyNetworkBridge.SYMBOLIC_NAME, String.class, BasicVisualLexicon.NODE_LABEL);
		vis.addVisualMappingFunction(nodeLabelMapping);	
		
		//edge default style
		
		//edge basic mappings
		//not showing edge labels by default for now
//		PassthroughMapping<String, ?> edgeLabelMapping = (PassthroughMapping<String, ?>)vmffP.createVisualMappingFunction(CyEdge.INTERACTION, String.class, BasicVisualLexicon.EDGE_LABEL);
//		vis.addVisualMappingFunction(edgeLabelMapping);
		
		VisualMappingFunctionFactory vmffD = services.getVisualMappingFunctionFactory("discrete");
		DiscreteMapping<Boolean, LineType> edgeLineMapping = (DiscreteMapping<Boolean, LineType>)vmffD.createVisualMappingFunction(CyNetworkBridge.NEGATED, Boolean.class, BasicVisualLexicon.EDGE_LINE_TYPE);
		edgeLineMapping.putMapValue(true, LineTypeVisualProperty.DASH_DOT);
		edgeLineMapping.putMapValue(false, LineTypeVisualProperty.SOLID);
		vis.addVisualMappingFunction(edgeLineMapping);
		
		ExtraDiffanyVisMappings extras = new ExtraDiffanyVisMappings(services);
		for (VisualMappingFunction<?, ?> mapping : extras.getMappings()){
			vis.addVisualMappingFunction(mapping);
		}
		
		
	}

	
	/**
	 * Add visual mappings according to the networks contained in the given {@link CyProject}
	 * 
	 * @param interactions list of edge types that should be mapped on color and shape
	 * @param ontology the {@link EdgeOntology} to determine visual properties of each edge
	 */
	public void updateInteractionMappings(Set<String> interactions, EdgeOntology ontology) {
		
		//get mapping factories
		VisualMappingFunctionFactory vmffD = services.getVisualMappingFunctionFactory("discrete");
		DiscreteMapping<String, Paint> edgeColorFunction = (DiscreteMapping<String, Paint>)vmffD.createVisualMappingFunction
				(CyEdge.INTERACTION, String.class, BasicVisualLexicon.EDGE_STROKE_UNSELECTED_PAINT);
		DiscreteMapping<String, Paint> edgeSelectedColorFunction = (DiscreteMapping<String, Paint>)vmffD.createVisualMappingFunction
				(CyEdge.INTERACTION, String.class, BasicVisualLexicon.EDGE_STROKE_SELECTED_PAINT);
		DiscreteMapping<String, ArrowShape> edgeTargetArrowFunction = (DiscreteMapping<String, ArrowShape>)vmffD.createVisualMappingFunction(CyEdge.INTERACTION, String.class, 
				BasicVisualLexicon.EDGE_TARGET_ARROW_SHAPE);
		
		//this is overwritten by each type of VisualStyle
		this.addInteractionMappings(interactions, ontology, 
				edgeColorFunction, 
				edgeSelectedColorFunction,
				edgeTargetArrowFunction);			

		
		vis.addVisualMappingFunction(edgeColorFunction);
		vis.addVisualMappingFunction(edgeSelectedColorFunction);
		vis.addVisualMappingFunction(edgeTargetArrowFunction);
		
	}
	
	/**
	 * Add mappings to the VizMapper according to the edges used in the project. All found edge interaction types
	 * are mapped to a visual style using the given {@link EdgeOntology}
	 * 
	 * @param interactionTypes list of edge types that should be mapped on color and shape
	 * @param edgeOntology the {@link EdgeOntology} to determine visual properties of each edge
	 * @param edgeColorFunction edge type to color mapping
	 * @param edgeSelectedColorFunction edge type to selected color mapping
	 * @param edgeTargetArrowFunction edge type to arrow head mapping
	 */
	protected abstract void addInteractionMappings(Set<String> interactionTypes, EdgeOntology edgeOntology, 
			DiscreteMapping<String, Paint> edgeColorFunction, 
			DiscreteMapping<String, Paint> edgeSelectedColorFunction,
			DiscreteMapping<String, ArrowShape> edgeTargetArrowFunction);

	/**
	 * 
	 * @return the actual {@link VisualStyle} object.
	 */
	public VisualStyle getVisualStyle(){
		return this.vis;
	}
	
	/**
	 * Convenience method to apply this VisualStyle to a CyNetworkView
	 * 
	 * @param view the CyNetworkView
	 */
	public void apply(CyNetworkView view){
		this.vis.apply(view);
		services.getVisualMappingManager().setVisualStyle(vis, view);
	}

}
