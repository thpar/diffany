package be.svlandeg.diffany.cytoscape.vizmapper;

import java.awt.Color;
import java.awt.Paint;
import java.util.Set;

import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.presentation.property.BasicVisualLexicon;
import org.cytoscape.view.presentation.property.NodeShapeVisualProperty;
import org.cytoscape.view.vizmap.VisualMappingFunctionFactory;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.view.vizmap.mappings.DiscreteMapping;
import org.cytoscape.view.vizmap.mappings.PassthroughMapping;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.CyProject.CyNetworkPair;
import be.svlandeg.diffany.internal.Services;
import be.svlandeg.diffany.semantics.EdgeOntology;

public abstract class AbstractVisualDiffanyStyle {

	protected String name;
	protected VisualStyle vis;
	protected Services services;
	
	//color definitions
	private static final Color NETWORK_BACKGROUND_COLOR = new Color(190, 195, 232);
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
		
		//node basic mappings
		VisualMappingFunctionFactory vmffP = services.getVisualMappingFunctionFactory("passthrough");
		PassthroughMapping<String, ?> nodeLabelMapping = (PassthroughMapping<String, ?>)vmffP.createVisualMappingFunction(CyNetwork.NAME, String.class, BasicVisualLexicon.NODE_LABEL);
		vis.addVisualMappingFunction(nodeLabelMapping);	
		
		//edge default style
		
		//edge basic mappings
		PassthroughMapping<String, ?> edgeLabelMapping = (PassthroughMapping<String, ?>)vmffP.createVisualMappingFunction(CyEdge.INTERACTION, String.class, BasicVisualLexicon.EDGE_LABEL);
		vis.addVisualMappingFunction(edgeLabelMapping);
	}

	
	/**
	 * Add visual mappings according to the networks contained in the given {@link CyProject}
	 */
	public void updateInteractionMappings(CyProject cyProject) {
		VisualMappingFunctionFactory vmffD = services.getVisualMappingFunctionFactory("discrete");
		DiscreteMapping<String, Paint> edgeColorFunction = (DiscreteMapping<String, Paint>)vmffD.createVisualMappingFunction
				(CyEdge.INTERACTION, String.class, BasicVisualLexicon.EDGE_STROKE_UNSELECTED_PAINT);
		DiscreteMapping<String, Paint> edgeSelectedColorFunction = (DiscreteMapping<String, Paint>)vmffD.createVisualMappingFunction
				(CyEdge.INTERACTION, String.class, BasicVisualLexicon.EDGE_STROKE_SELECTED_PAINT);
		
		EdgeOntology edgeOntology = cyProject.getEdgeOntology();
				
		CyNetwork refNet = cyProject.getReferenceNetwork();
		Set<CyNetwork> conditionNetworks = cyProject.getConditionalNetworks();
		Set<CyNetworkPair> resultNets = cyProject.getResultNetworks();
		
		this.addInteractionMappings(refNet, edgeOntology, edgeColorFunction, edgeSelectedColorFunction);
		
		for (CyNetwork condNet : conditionNetworks){
			this.addInteractionMappings(condNet, edgeOntology, edgeColorFunction, edgeSelectedColorFunction);			
		}
		
		for (CyNetworkPair pair : resultNets){
			this.addInteractionMappings(pair.diffNet, edgeOntology, edgeColorFunction, edgeSelectedColorFunction);
			this.addInteractionMappings(pair.overlapNet, edgeOntology, edgeColorFunction, edgeSelectedColorFunction);
		}

		
		vis.addVisualMappingFunction(edgeColorFunction);
		vis.addVisualMappingFunction(edgeSelectedColorFunction);
		
	}
	
	/**
	 * Add mappings to the VizMapper according to the edges used in the project. All found edge interaction types
	 * are mapped to a visual style using the given {@link EdgeOntology}
	 * 
	 * @param refNet
	 * @param edgeOntology
	 * @param edgeColorFunction
	 */
	protected abstract void addInteractionMappings(CyNetwork refNet, EdgeOntology edgeOntology, 
			DiscreteMapping<String, Paint> edgeColorFunction, DiscreteMapping<String, Paint> edgeSelectedColorFunction);

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
	}

}
