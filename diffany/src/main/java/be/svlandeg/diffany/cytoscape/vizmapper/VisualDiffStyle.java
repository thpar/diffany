package be.svlandeg.diffany.cytoscape.vizmapper;

import java.awt.Paint;
import java.util.Set;

import org.cytoscape.model.CyEdge;
import org.cytoscape.view.vizmap.mappings.DiscreteMapping;

import be.svlandeg.diffany.internal.Services;
import be.svlandeg.diffany.semantics.EdgeOntology;

public class VisualDiffStyle extends AbstractVisualDiffanyStyle {

	public VisualDiffStyle(Services services) {
		super("Diffany - Differential", services);
	}

	
	@Override
	protected void addInteractionMappings(Set<String> interactionTypes, EdgeOntology edgeOntology, 
			DiscreteMapping<String, Paint> edgeColorFunction, DiscreteMapping<String, Paint> edgeSelectedColorFunction) {
		System.out.println("Updating Differential Visual Style");
		for (String type : interactionTypes) {
			Paint paint = edgeOntology.getDifferentialEdgeStyle(type);
			edgeColorFunction.putMapValue(type, paint);
			System.out.println("Mapped "+type+" to "+paint);
		}
	}
	

}
