package be.svlandeg.diffany.cytoscape.gui;

import java.util.ArrayList;
import java.util.List;
import java.util.Observable;
import java.util.Observer;

import javax.swing.table.AbstractTableModel;

import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CySubNetwork;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.NetworkEntry;

/**
 * The model for the list of available networks in the selected collection.
 * 
 * @author thpar
 *
 */
public class SelectionTableModel extends AbstractTableModel{

	private static final long serialVersionUID = 1L;

	
	/**
	 * Networks to be listed in the selection list
	 */
	private List<NetworkEntry> networkEntries = new ArrayList<NetworkEntry>();
	
	private int referenceRow = 0;

	/**
	 * Column headers
	 */
	private String[] columns = {"Include", "Network", "Reference"};

	
	
	/**
	 * Create a new model based on the general {@link Model} and add it as an {@link Observer} the contained
	 * {@link GUIModel}.
	 * @param model
	 */
	public SelectionTableModel() {
	}
	
	
	@Override
	public int getRowCount() {
		return networkEntries.size();
	}

	@Override
	public int getColumnCount() {
		return columns.length;
	}

	@Override
	public String getColumnName(int columnIndex) {
		return columns[columnIndex];
	}

	@Override
	public Class<?> getColumnClass(int columnIndex) {
		switch(columnIndex){
		case 0:
		case 2:
			return Boolean.class;
		case 1: 
		default:
			return String.class;
		}
	}

	@Override
	public boolean isCellEditable(int rowIndex, int columnIndex) {
		NetworkEntry entry = networkEntries.get(rowIndex);
		switch(columnIndex){
		case 1:
			return false;
		case 0:
			return !entry.isReference();
		case 2:
			return entry.isSelected() && !entry.isReference();
		}
		return false;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		NetworkEntry entry = networkEntries.get(rowIndex);
		switch(columnIndex){
		case 0:
			return entry.isSelected();
		case 1:
			return entry.getName();
		case 2:
			return entry.isReference();
		}
		return false;
	}

	@Override
	public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
		NetworkEntry entry = networkEntries.get(rowIndex);
		boolean check = (Boolean)aValue;
		switch(columnIndex){
		case 0:
			entry.setSelected(check);
			break;
		case 2:
			this.setReference(rowIndex);
			break;
		}
	}

	private void setReference(int row){
		int oldRefRow = this.referenceRow;
		NetworkEntry oldRefEntry = this.networkEntries.get(oldRefRow);
		oldRefEntry.setReference(false);
		
		NetworkEntry newRefEntry = this.networkEntries.get(row);
		newRefEntry.setReference(true);
		this.referenceRow = row;
		
		this.fireTableRowsUpdated(oldRefRow, oldRefRow);
		
	}
	
	/**
	 * Reload the {@link NetworkEntry}s based on the subnetworks from the selected network collections.
	 * 
	 * @param list of {@link CySubNetwork}s to be displayed.
	 */
	public void refresh(List<CySubNetwork> subNets) {
		networkEntries = new ArrayList<NetworkEntry>();
				
		for (CySubNetwork subNet : subNets){
			NetworkEntry entry = new NetworkEntry(subNet);
			entry.setSelected(true);
			entry.setReference(false);
			networkEntries.add(entry);
		}

		if (networkEntries.size() > 0){
			this.setReference(0);
		}
		this.fireTableDataChanged();
	}
	
	
}
