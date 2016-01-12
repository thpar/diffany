package be.svlandeg.diffany.cytoscape.gui;

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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import javax.swing.table.AbstractTableModel;

import be.svlandeg.diffany.cytoscape.CyProject;

/**
 * Model that controls the content of the edge (interaction) filter
 * 
 * @author Thomas Van Parys
 *
 */
public class EdgeFilterTableModel extends AbstractTableModel {


	private static final long serialVersionUID = 1L;
	private String[] columns = {"Interaction", "Display"};
	
	private Map<String, Boolean> sourceInteractions = new HashMap<String,Boolean>();
	private Map<String, Boolean> diffInteractions = new HashMap<String,Boolean>();
	
	
	@Override
	public int getRowCount() {
		int total = sourceInteractions.size();
		if (diffInteractions.size() !=0){
			total = total + diffInteractions.size();
		}
		return total;
	}

	@Override
	public int getColumnCount() {
		return columns.length;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		switch(columnIndex){
		case 0:
			return getInteractionInRow(rowIndex);
		case 1:
			return getValueInRow(rowIndex);
		default:
			return null;
		}
	}
	
	
	/**
	 * The interactions are sorted alphabetically in the table. Get the one 
	 * with given row index.
	 * 
	 * @param rowIndex
	 * @return interaction name at given row index
	 */
	private String getInteractionInRow(int rowIndex){
		if (rowIndex<sourceInteractions.size()){
			//return value from Source list
			List<String> interactionNames = new ArrayList<String>(sourceInteractions.keySet());
			Collections.sort(interactionNames);			
			return interactionNames.get(rowIndex);
		} else{
			//return value from Diff list
			List<String> interactionNames = new ArrayList<String>(diffInteractions.keySet());
			Collections.sort(interactionNames);			
			return interactionNames.get(rowIndex-sourceInteractions.size());
		} 
	}
	/**
	 * The interactions are sorted alphabetically in the table. Get value of the one 
	 * with given row index.
	 * 
	 * @param rowIndex
	 * @return interaction name at given row index
	 */
	private boolean getValueInRow(int rowIndex){
		if (rowIndex<sourceInteractions.size()){
			return sourceInteractions.get(getInteractionInRow(rowIndex));
		} else {
			return diffInteractions.get(getInteractionInRow(rowIndex));
		}
	}
	
	@Override
	public String getColumnName(int columnIndex) {
		return columns[columnIndex];
	}
	
	@Override
	public Class<?> getColumnClass(int columnIndex) {
		switch(columnIndex){
		case 1:
			return Boolean.class;
		case 0:
		default:
			return String.class;
		}
	}
	
	@Override
	public boolean isCellEditable(int rowIndex, int columnIndex) {
		switch(columnIndex){
		case 1:
			return true;
		case 0:
		default:
			return false;
		}
		
	}
	
	@Override
	public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
		if (columnIndex == 1){
			boolean value = (Boolean)aValue;
			String interaction = this.getInteractionInRow(rowIndex);
			if (rowIndex<sourceInteractions.size()){
				boolean currentValue = sourceInteractions.get(interaction);
				if (currentValue != value){
					sourceInteractions.put(interaction, value);
					this.fireTableDataChanged();
				}				
			} else {
				boolean currentValue = diffInteractions.get(interaction);
				if (currentValue != value){
					diffInteractions.put(interaction, value);
					this.fireTableDataChanged();
				}				
			}
		}
	}

	
	/**
	 * Repopulate the model with all interactions found in the currently selected {@link CyProject} 
	 * 
	 * @param selectedProject the currently selected project
	 */
	public void refresh(CyProject selectedProject) {
		Set<String> projectSourceInteractions = selectedProject.getAllSourceInteractions();		
		Map<String, Boolean> newSourceInteractions = new HashMap<String, Boolean>();
		for (String interaction : projectSourceInteractions){
			boolean value = true;
			if (this.sourceInteractions.containsKey(interaction)){
				value = this.sourceInteractions.get(interaction);
			}
			newSourceInteractions.put(interaction, value);
		}
		this.sourceInteractions = newSourceInteractions;
		
		Set<String> projectDiffInteractions = selectedProject.getAllDifferentialInteractions();
		Map<String, Boolean> newDiffInteractions = new HashMap<String, Boolean>();
		for (String interaction : projectDiffInteractions){
			boolean value = true;
			if (this.diffInteractions.containsKey(interaction)){
				value = this.diffInteractions.get(interaction);
			}
			newDiffInteractions.put(interaction, value);
		}
		this.diffInteractions = newDiffInteractions;
		
		this.fireTableDataChanged();
	}
	
	/**
	 * Get a set of edge types (interactions) that should be hidden.
	 * @return the set of hidden edge types 
	 */
	public Set<String> getHiddenInteractions(){
		Set<String> hidden = new HashSet<String>();
		for (Entry<String, Boolean> interaction : sourceInteractions.entrySet()){
			if (!interaction.getValue()){
				hidden.add(interaction.getKey());
			}
		}
		for (Entry<String, Boolean> interaction : diffInteractions.entrySet()){
			if (!interaction.getValue()){
				hidden.add(interaction.getKey());
			}
		}
		return hidden;
	}

	/**
	 * Empty the table and the complete table model
	 */
	public void clear() {
		this.sourceInteractions = new HashMap<String, Boolean>();
		this.diffInteractions = new HashMap<String, Boolean>();
		this.fireTableDataChanged();
	}

}
