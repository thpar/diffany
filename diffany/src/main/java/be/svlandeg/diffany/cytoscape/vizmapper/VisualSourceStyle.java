package be.svlandeg.diffany.cytoscape.vizmapper;

import java.awt.Paint;

import org.cytoscape.view.presentation.property.BasicVisualLexicon;
import org.cytoscape.view.vizmap.VisualMappingFunction;
import org.cytoscape.view.vizmap.VisualMappingFunctionFactory;
import org.cytoscape.view.vizmap.mappings.PassthroughMapping;

import be.svlandeg.diffany.internal.Services;

public class VisualSourceStyle extends VisualDiffanyStyle {
	
	static final String INTERACTION = "interaction";
	
	
	public VisualSourceStyle(Services services) {
		super("Diffany - Source", services);
	}

	@Override
	protected void init() {
		VisualMappingFunctionFactory vmffP = services.getVisualMappingFunctionFactory("passthrough");
		
		PassthroughMapping<String, ?> mapping = (PassthroughMapping<String, ?>)vmffP.createVisualMappingFunction("interaction", String.class, BasicVisualLexicon.EDGE_LABEL);
		
		vis.addVisualMappingFunction(mapping);
		
		
		VisualMappingFunctionFactory vmffD = services.getVisualMappingFunctionFactory("discrete");
		
		VisualMappingFunction<String, Paint> edgeColorFunction = vmffD.createVisualMappingFunction(INTERACTION, String.class, BasicVisualLexicon.EDGE_PAINT);
		
	}

}
