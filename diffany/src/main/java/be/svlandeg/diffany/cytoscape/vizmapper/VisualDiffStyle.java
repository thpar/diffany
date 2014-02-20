package be.svlandeg.diffany.cytoscape.vizmapper;

import java.awt.Color;
import java.awt.Paint;
import java.util.Set;

import org.cytoscape.view.presentation.property.ArrowShapeVisualProperty;
import org.cytoscape.view.presentation.property.BasicVisualLexicon;
import org.cytoscape.view.presentation.property.values.ArrowShape;
import org.cytoscape.view.vizmap.VisualMappingFunctionFactory;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.view.vizmap.mappings.BoundaryRangeValues;
import org.cytoscape.view.vizmap.mappings.ContinuousMapping;
import org.cytoscape.view.vizmap.mappings.DiscreteMapping;

import be.svlandeg.diffany.cytoscape.CyNetworkBridge;
import be.svlandeg.diffany.cytoscape.internal.Services;
import be.svlandeg.diffany.semantics.EdgeOntology;
import be.svlandeg.diffany.visualstyle.EdgeDrawing;
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
		
		EdgeDrawing diffDrawing = edgeOntology.getDifferentialEdgeDrawing();
		
		//map weight to edge width
		VisualMappingFunctionFactory vmffC = services.getVisualMappingFunctionFactory("continuous");
		ContinuousMapping<Double, Double> edgeWidthMapping = (ContinuousMapping<Double, Double>)vmffC.createVisualMappingFunction(CyNetworkBridge.WEIGHT, Double.class, BasicVisualLexicon.EDGE_WIDTH);
		edgeWidthMapping.addPoint(diffDrawing.getMinWeight(), new BoundaryRangeValues<Double>(0.0d,0.0d,0.0d));
		edgeWidthMapping.addPoint(diffDrawing.getMaxWeight(), new BoundaryRangeValues<Double>(20d,20d,20d));
		vis.addVisualMappingFunction(edgeWidthMapping);
		
		for (String type : interactionTypes) {
			EdgeStyle edgeStyle = diffDrawing.getEdgeStyle(type);
		
			
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
