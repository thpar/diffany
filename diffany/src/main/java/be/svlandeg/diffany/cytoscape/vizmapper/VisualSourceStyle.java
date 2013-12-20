package be.svlandeg.diffany.cytoscape.vizmapper;

import java.awt.Color;
import java.awt.Paint;

import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.vizmap.mappings.DiscreteMapping;

import be.svlandeg.diffany.internal.Services;
import be.svlandeg.diffany.semantics.EdgeOntology;

public class VisualSourceStyle extends AbstractVisualDiffanyStyle {
	
	
	
	public VisualSourceStyle(Services services) {
		super("Diffany - Source", services);
	}

	@Override
	protected void addInteractionMappings(CyNetwork cyNetwork, EdgeOntology edgeOntology, 
			DiscreteMapping<String, Paint> edgeColorFunction,
			DiscreteMapping<String, Paint> edgeSelectedColorFunction) {
		System.out.println("Updating Source Visual Style");
		for (CyEdge cyEdge : cyNetwork.getEdgeList()) {
			String interaction = cyNetwork.getRow(cyEdge).get(cyEdge.INTERACTION, String.class);
			
			Color paint = (Color)edgeOntology.getSourceEdgeStyle(interaction);
			edgeColorFunction.putMapValue(interaction, paint);
			
			Color darkPaint = paint.darker().darker();
			edgeSelectedColorFunction.putMapValue(interaction, darkPaint);
			System.out.println("Mapped "+interaction+" to "+paint);
		}
	}

	

}
