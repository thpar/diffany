package be.svlandeg.diffany.cytoscape.vizmapper;

import java.awt.Color;
import java.awt.Paint;
import java.util.Set;

import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.presentation.property.BasicVisualLexicon;
import org.cytoscape.view.presentation.property.LineTypeVisualProperty;
import org.cytoscape.view.presentation.property.NodeShapeVisualProperty;
import org.cytoscape.view.presentation.property.values.ArrowShape;
import org.cytoscape.view.presentation.property.values.LineType;
import org.cytoscape.view.vizmap.VisualMappingFunctionFactory;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.view.vizmap.mappings.BoundaryRangeValues;
import org.cytoscape.view.vizmap.mappings.ContinuousMapping;
import org.cytoscape.view.vizmap.mappings.DiscreteMapping;
import org.cytoscape.view.vizmap.mappings.PassthroughMapping;

import be.svlandeg.diffany.cytoscape.CyNetworkBridge;
import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.internal.Services;
import be.svlandeg.diffany.semantics.EdgeOntology;

/**
 * A {@link VisualStyle} wrapper with default values for all Diffany styles. All new styles within this app
 * should extend this class.
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
	private static final Color NODE_COLOR = Color.YELLOW;
	
	
	
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
		vis.setDefaultValue(BasicVisualLexicon.NODE_SHAPE, NodeShapeVisualProperty.ROUND_RECTANGLE);
		vis.setDefaultValue(BasicVisualLexicon.NODE_FILL_COLOR, NODE_COLOR);
		vis.setDefaultValue(BasicVisualLexicon.NODE_SIZE, 30d);
		
		//node basic mappings
		VisualMappingFunctionFactory vmffP = services.getVisualMappingFunctionFactory("passthrough");
		PassthroughMapping<String, ?> nodeLabelMapping = (PassthroughMapping<String, ?>)vmffP.createVisualMappingFunction(CyNetwork.NAME, String.class, BasicVisualLexicon.NODE_LABEL);
		vis.addVisualMappingFunction(nodeLabelMapping);	
		
		//edge default style
		
		//edge basic mappings
		PassthroughMapping<String, ?> edgeLabelMapping = (PassthroughMapping<String, ?>)vmffP.createVisualMappingFunction(CyEdge.INTERACTION, String.class, BasicVisualLexicon.EDGE_LABEL);
		vis.addVisualMappingFunction(edgeLabelMapping);
		
		VisualMappingFunctionFactory vmffD = services.getVisualMappingFunctionFactory("discrete");
		DiscreteMapping<Boolean, LineType> edgeLineMapping = (DiscreteMapping<Boolean, LineType>)vmffD.createVisualMappingFunction(CyNetworkBridge.NEGATED, Boolean.class, BasicVisualLexicon.EDGE_LINE_TYPE);
		edgeLineMapping.putMapValue(true, LineTypeVisualProperty.DASH_DOT);
		edgeLineMapping.putMapValue(false, LineTypeVisualProperty.SOLID);
		vis.addVisualMappingFunction(edgeLineMapping);		
		
		
	}

	
	/**
	 * Add visual mappings according to the networks contained in the given {@link CyProject}
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
	 * @param interactionTypes All interactions used in the networks of the project.
	 * @param edgeOntology
	 * @param edgeColorFunction
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
	 * Convenience method to apply this VisualStyle to a CyView
	 * 
	 * @param view
	 */
	public void apply(CyNetworkView view){
		this.vis.apply(view);
		services.getVisualMappingManager().setVisualStyle(vis, view);
	}

}
