package be.svlandeg.diffany.cytoscape.gui;

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

public class EdgeFilterTableModel extends AbstractTableModel {


	private static final long serialVersionUID = 1L;
	private String[] columns = {"Interaction", "Display"};
	
	private Map<String, Boolean> sourceInteractions = new HashMap<String,Boolean>();
	private Map<String, Boolean> diffInteractions = new HashMap<String,Boolean>();
	
	
	@Override
	public int getRowCount() {
		int total = sourceInteractions.size();
		if (diffInteractions.size() !=0){
			total = total + diffInteractions.size()+1;
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
		} else if (rowIndex>sourceInteractions.size()){
			//return value from Diff list
			List<String> interactionNames = new ArrayList<String>(diffInteractions.keySet());
			Collections.sort(interactionNames);			
			return interactionNames.get(rowIndex-sourceInteractions.size()-1);
		} else {
			//return separator row
			return new String();
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
		} else if (rowIndex>sourceInteractions.size()){
			return diffInteractions.get(getInteractionInRow(rowIndex));
		} else {
			//return separator row
			return false;
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
		if (rowIndex==sourceInteractions.size()){
			return false;
		}
		
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
			} else if (rowIndex > sourceInteractions.size()){
				boolean currentValue = diffInteractions.get(interaction);
				if (currentValue != value){
					diffInteractions.put(interaction, value);
					this.fireTableDataChanged();
				}				
			}
		}
	}

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

	public void clear() {
		this.sourceInteractions = new HashMap<String, Boolean>();
		this.diffInteractions = new HashMap<String, Boolean>();
		this.fireTableDataChanged();
	}

}
