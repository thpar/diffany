package be.svlandeg.diffany.cytoscape.vizmapper;

import java.awt.Color;
import java.awt.Paint;
import java.util.Set;

import org.cytoscape.view.presentation.property.ArrowShapeVisualProperty;
import org.cytoscape.view.presentation.property.values.ArrowShape;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.view.vizmap.mappings.DiscreteMapping;

import be.svlandeg.diffany.internal.Services;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.visualstyle.EdgeStyle;
import be.svlandeg.diffany.visualstyle.EdgeStyle.ArrowHead;

/**
 * {@link VisualStyle} to be applied on differential networks.
 * 
 * @author Thomas Van Parys
 *
 */
public class VisualDiffStyle extends AbstractVisualDiffanyStyle {

	/**
	 * Construct and register new visual style for differential networks.
	 * 
	 * @param services app services
	 */
	public VisualDiffStyle(Services services) {
		super("Diffany - Differential", services);
	}

	
	@Override
	protected void addInteractionMappings(Set<String> interactionTypes, EdgeOntology edgeOntology, 
			DiscreteMapping<String, Paint> edgeColorFunction, 
			DiscreteMapping<String, Paint> edgeSelectedColorFunction,
			DiscreteMapping<String, ArrowShape> edgeTargetArrowFunction) {
		for (String type : interactionTypes) {
			EdgeStyle edgeStyle = edgeOntology.getDifferentialEdgeDrawing().getEdgeStyle(type);
			
			//edge color
			Color paint = edgeStyle.getColor();
			edgeColorFunction.putMapValue(type, paint);
			Color darkPaint = paint.darker().darker();
			edgeSelectedColorFunction.putMapValue(type, darkPaint);
			
			//arrow head
			ArrowHead arrowHead = edgeStyle.getArrowHead();
			ArrowShape shape;
			switch(arrowHead){
			default:
			case NONE:
				shape = ArrowShapeVisualProperty.NONE;
				break;
			case ARROW:
				shape = ArrowShapeVisualProperty.ARROW;
				break;
			case T:
				shape = ArrowShapeVisualProperty.T;
				break;
			case DIAMOND:
				shape = ArrowShapeVisualProperty.DIAMOND;
				break;
			}
			edgeTargetArrowFunction.putMapValue(type,shape);
			
			//edge stroke
			
			
		}
	}
	

}
