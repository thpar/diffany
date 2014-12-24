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
	 * @return the extra visual mappings
	 */
	public List<VisualMappingFunction<?,?>> getMappings();
	
}
