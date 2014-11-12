package be.svlandeg.diffany.cytoscape.vizmapper;

import java.awt.Color;
import java.awt.Paint;
import java.util.ArrayList;
import java.util.List;

import org.cytoscape.view.presentation.property.BasicVisualLexicon;
import org.cytoscape.view.presentation.property.NodeShapeVisualProperty;
import org.cytoscape.view.presentation.property.values.NodeShape;
import org.cytoscape.view.vizmap.VisualMappingFunction;
import org.cytoscape.view.vizmap.VisualMappingFunctionFactory;
import org.cytoscape.view.vizmap.mappings.DiscreteMapping;

import be.svlandeg.diffany.cytoscape.internal.Services;



/**
 * Define extra visual mappings that are not specifically needed by the
 * Diffany algorithm.
 * 
 * @author thpar
 *
 */
public class ExtraDiffanyVisMappings implements ExtraVisMappings{
	List<VisualMappingFunction<?,?>> mappings = new ArrayList<VisualMappingFunction<?,?>>();

	private final String PHOSPHORYLATION_COLUMN = "phosphorylation_site";
	private final NodeShape PHOSPHORYLATION_SHAPE = NodeShapeVisualProperty.DIAMOND;
	private final String PHOSPHORYLATION_YES = "yes";
	private final String PHOSPHORYLATION_NO = "no";
	
	private final String KINASE_COLUMN = "kinase_function";
	private final double KINASE_BORDER_SIZE = 2d;
	private final Color KINASE_BORDER_COLOR = Color.BLACK;
	private final String KINASE_YES = "yes";
	private final String KINASE_NO = "no";
	
	private final String DE_COLUMN = "differentially_expressed";	
	private final Color DE_UP_COLOR = new Color(255,255,15);      //yellow
	private final Color DE_DOWN_COLOR = new Color(62,181,255);    //blue
	private final Color DE_NEUTRAL_COLOR = Color.GRAY; //gray
	private final String DE_UP_VALUE = "up-regulated"; 
	private final String DE_DOWN_VALUE = "down-regulated"; 
	private final String DE_NEUTRAL_VALUE = "no"; 
	
	
	public ExtraDiffanyVisMappings(Services services) {
		VisualMappingFunctionFactory vmffD = services.getVisualMappingFunctionFactory("discrete");
		DiscreteMapping<String, NodeShape> nodeShapeMapping = 
				(DiscreteMapping<String, NodeShape>)vmffD.createVisualMappingFunction(PHOSPHORYLATION_COLUMN, String.class, BasicVisualLexicon.NODE_SHAPE);
		nodeShapeMapping.putMapValue(PHOSPHORYLATION_YES, PHOSPHORYLATION_SHAPE);
		mappings.add(nodeShapeMapping);
		
		DiscreteMapping<String, Double> nodeBorderMapping = 
				(DiscreteMapping<String, Double>)vmffD.createVisualMappingFunction(KINASE_COLUMN, String.class, BasicVisualLexicon.NODE_BORDER_WIDTH);
		nodeBorderMapping.putMapValue(KINASE_YES, KINASE_BORDER_SIZE);
		mappings.add(nodeBorderMapping);
		
		DiscreteMapping<String, Paint> nodeBorderColorMapping = 
				(DiscreteMapping<String, Paint>)vmffD.createVisualMappingFunction(KINASE_COLUMN, String.class, BasicVisualLexicon.NODE_BORDER_PAINT);
		nodeBorderColorMapping.putMapValue(KINASE_YES, KINASE_BORDER_COLOR);
		mappings.add(nodeBorderColorMapping);
		
		DiscreteMapping<String, Paint> nodeColorMapping = 
				(DiscreteMapping<String, Paint>)vmffD.createVisualMappingFunction(DE_COLUMN, String.class, BasicVisualLexicon.NODE_FILL_COLOR);
		nodeBorderMapping.putMapValue(KINASE_YES, KINASE_BORDER_SIZE);
		nodeColorMapping.putMapValue(DE_UP_VALUE, DE_UP_COLOR);
		nodeColorMapping.putMapValue(DE_DOWN_VALUE, DE_DOWN_COLOR);
		mappings.add(nodeColorMapping);
		
	}
	
	@Override
	public List<VisualMappingFunction<?,?>> getMappings() {
		return mappings;
	}

}
