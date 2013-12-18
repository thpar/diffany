package be.svlandeg.diffany.cytoscape.vizmapper;

import java.awt.Paint;

import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.vizmap.mappings.DiscreteMapping;

import be.svlandeg.diffany.internal.Services;
import be.svlandeg.diffany.semantics.EdgeOntology;

public class VisualDiffStyle extends AbstractVisualDiffanyStyle {

	public VisualDiffStyle(Services services) {
		super("Diffany - Differential", services);
	}

	
	@Override
	protected void addInteractionMappings(CyNetwork cyNetwork, EdgeOntology edgeOntology, DiscreteMapping<String, Paint> edgeColorFunction) {
		for (CyEdge cyEdge : cyNetwork.getEdgeList()) {
			String interaction = cyNetwork.getRow(cyEdge).get(cyEdge.INTERACTION, String.class);
			Paint paint = edgeOntology.getDifferentialEdgeStyle(interaction);
			edgeColorFunction.putMapValue(interaction, paint);
		}
	}
	

}
