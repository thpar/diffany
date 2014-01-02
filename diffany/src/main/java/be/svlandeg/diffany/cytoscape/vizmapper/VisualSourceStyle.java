package be.svlandeg.diffany.cytoscape.vizmapper;

import java.awt.Color;
import java.awt.Paint;
import java.util.Set;

import org.cytoscape.model.CyEdge;
import org.cytoscape.view.vizmap.mappings.DiscreteMapping;

import be.svlandeg.diffany.internal.Services;
import be.svlandeg.diffany.semantics.EdgeOntology;

public class VisualSourceStyle extends AbstractVisualDiffanyStyle {
	
	
	
	public VisualSourceStyle(Services services) {
		super("Diffany - Source", services);
	}

	@Override
	protected void addInteractionMappings(Set<String> interactionTypes, EdgeOntology edgeOntology, 
			DiscreteMapping<String, Paint> edgeColorFunction,
			DiscreteMapping<String, Paint> edgeSelectedColorFunction) {
		
		System.out.println("Updating Source Visual Style");
		for (String type : interactionTypes) {
			
			Color paint = (Color)edgeOntology.getSourceEdgeStyle(type);
			edgeColorFunction.putMapValue(type, paint);
			
			Color darkPaint = paint.darker().darker();
			edgeSelectedColorFunction.putMapValue(type, darkPaint);
			System.out.println("Mapped "+type+" to "+paint);
		}
	}

	

}
