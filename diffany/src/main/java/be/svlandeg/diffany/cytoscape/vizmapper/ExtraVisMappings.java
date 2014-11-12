package be.svlandeg.diffany.cytoscape.vizmapper;

import java.util.List;

import org.cytoscape.view.vizmap.VisualMappingFunction;

/**
 * Define extra visual mappings that are not specifically needed by the
 * Diffany algorithm.
 * @author thpar
 *
 */
public interface ExtraVisMappings {

	/**
	 * Get all the extra visual mappings. These mappings will be loaded at the end of every VizMapper refresh,
	 * so they will be overriding any previous defined mappings that might have been defined.
	 * 
	 * 
	 * @return
	 */
	public List<VisualMappingFunction<?,?>> getMappings();
	
}
