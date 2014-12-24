package be.svlandeg.diffany.cytoscape.vizmapper;

/*
 * #%L
 * Diffany
 * %%
 * Copyright (C) 2014 PSB/UGent - Sofie Van Landeghem and Thomas Van Parys
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */

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
	
	/**
	 * Define the extra VizMapper settings
	 * 
	 * @param services App services to supply the mapping factories
	 */
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
		nodeColorMapping.putMapValue(DE_UP_VALUE, DE_UP_COLOR);
		nodeColorMapping.putMapValue(DE_DOWN_VALUE, DE_DOWN_COLOR);
		mappings.add(nodeColorMapping);
		
		DiscreteMapping<String, Paint> nodeSelectedColorMapping = 
				(DiscreteMapping<String, Paint>)vmffD.createVisualMappingFunction(DE_COLUMN, String.class, BasicVisualLexicon.NODE_SELECTED_PAINT);
		nodeSelectedColorMapping.putMapValue(DE_UP_VALUE, DE_UP_COLOR.darker().darker());
		nodeSelectedColorMapping.putMapValue(DE_DOWN_VALUE, DE_DOWN_COLOR.darker().darker());
		mappings.add(nodeSelectedColorMapping);
		
	}
	
	@Override
	public List<VisualMappingFunction<?,?>> getMappings() {
		return mappings;
	}

}
