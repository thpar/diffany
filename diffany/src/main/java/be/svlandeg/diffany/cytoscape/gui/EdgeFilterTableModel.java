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
	
	private Map<String, Boolean> interactions = new HashMap<String,Boolean>();
		
	
	@Override
	public int getRowCount() {
		return interactions.size();
	}

	@Override
	public int getColumnCount() {
		return columns.length;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		String name = getInteractionInRow(rowIndex);
		switch(columnIndex){
		case 0:
			return name;
		case 1:
			return interactions.get(name);
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
		List<String> interactionNames = new ArrayList<String>(interactions.keySet());
		Collections.sort(interactionNames);
		return interactionNames.get(rowIndex);
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
			boolean currentValue = interactions.get(interaction);
			if (currentValue != value){
				interactions.put(interaction, value);
				this.fireTableDataChanged();
			}
		}
	}

	public void refresh(CyProject selectedProject) {
		Set<String> projectInteractions = selectedProject.getAllInteractions();
		
		Map<String, Boolean> newInteractions = new HashMap<String, Boolean>();
		for (String interaction : projectInteractions){
			boolean value = true;
			if (this.interactions.containsKey(interaction)){
				value = this.interactions.get(interaction);
			}
			newInteractions.put(interaction, value);
		}
		this.interactions = newInteractions;
		
		System.out.println("---");
		for (Entry<String, Boolean> entry : this.interactions.entrySet()){
			System.out.println(entry.getKey()+" = "+entry.getValue());
		}
		System.out.println("---");
		
		this.fireTableDataChanged();
	}
	
	/**
	 * Get a set of edge types (interactions) that should be hidden.
	 */
	public Set<String> getHiddenInteractions(){
		Set<String> hidden = new HashSet<String>();
		for (Entry<String, Boolean> interaction : interactions.entrySet()){
			if (!interaction.getValue()){
				hidden.add(interaction.getKey());
			}
		}
		return hidden;
	}

	public void clear() {
		this.interactions = new HashMap<String, Boolean>();
		this.fireTableDataChanged();
	}

}
