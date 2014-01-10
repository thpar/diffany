package be.svlandeg.diffany.cytoscape.gui;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.table.AbstractTableModel;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyRow;
import org.cytoscape.model.subnetwork.CySubNetwork;

import be.svlandeg.diffany.cytoscape.NetworkEntry;

/**
 * The model for the list of available networks in the selected collection.
 * 
 * @author Thomas Van Parys
 *
 */
public class SelectionTableModel extends AbstractTableModel{

	private static final long serialVersionUID = 1L;

	
	/**
	 * Networks to be listed in the selection list
	 */
	private List<NetworkEntry> networkEntries = new ArrayList<NetworkEntry>();
	private Map<CyNetwork, NetworkEntry> entryMap = new HashMap<CyNetwork, NetworkEntry>();
	
	private int referenceRow = 0;

	/**
	 * Column headers
	 */
	private String[] columns = {"Include", "Network", "Reference"};


	
	/**
	 * Create a new table model.
	 * 
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
			return !entry.isReference();
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
			this.fireTableDataChanged();
			break;
		case 2:
			if (!entry.isSelected()){
				entry.setSelected(true);
			}
			this.setReference(rowIndex);
			break;
		}
	}

	private void setReference(int row){
		int oldRefRow = this.referenceRow;
		if (oldRefRow != -1){
			NetworkEntry oldRefEntry = this.networkEntries.get(oldRefRow);
			oldRefEntry.setReference(false);			
		}
		
		NetworkEntry newRefEntry = this.networkEntries.get(row);
		newRefEntry.setReference(true);
		this.referenceRow = row;
		
		this.fireTableDataChanged();
		
	}
	

	/**
	 * Reload the {@link NetworkEntry}s based on the subnetworks from the selected network collections.
	 * Compare with the previous list to maintain selections.
	 * 
	 * @param list of {@link CySubNetwork}s to be displayed.
	 */
	public void refresh(List<CySubNetwork> subNets) {
		Map<CyNetwork, NetworkEntry> oldMap = entryMap;
		
		networkEntries = new ArrayList<NetworkEntry>();
		entryMap = new HashMap<CyNetwork, NetworkEntry>();		
		
		this.referenceRow = -1;
		for (CySubNetwork subNet : subNets){
			if (!isGhostNetwork(subNet)){
				if (oldMap.containsKey(subNet)){
					NetworkEntry origEntry = oldMap.get(subNet);
					networkEntries.add(origEntry);
					entryMap.put(subNet, origEntry);
					if (origEntry.isReference()){
						this.referenceRow = networkEntries.indexOf(origEntry);
					}
				} else {
					NetworkEntry entry = new NetworkEntry(subNet);
					entry.setSelected(false);
					entry.setReference(false);
					networkEntries.add(entry);
					entryMap.put(subNet, entry);
				}
				
			}
		}
		this.fireTableDataChanged();
	}
	
	/**
	 * Allows to ignore the ghost network that is often created when importing a network or loading a session.
	 * 
	 * @param subNet
	 * @return
	 */
	private boolean isGhostNetwork(CySubNetwork subNet) {
		String name = subNet.getRow(subNet).get(CyNetwork.NAME, String.class);
		return name==null || name.isEmpty();			
	}


	/**
	 * Get all conditional {@link CyNetwork}s in this project.
	 * 
	 * @return all conditional {@link CyNetwork}s in this project.
	 */
	public Set<CyNetwork> getConditionalNetworks(){
		Set<CyNetwork> condSet = new HashSet<CyNetwork>();
		for (NetworkEntry network : networkEntries){
			if (network.isSelected() && !network.isReference()){
				condSet.add(network.getNetwork());
			}
		}
		return condSet;
	}
	
	/**
	 * Gets the {@link CyNetwork} to be used as reference network.
	 * 
	 * @return the {@link CyNetwork} to be used as reference network. Returns null is no reference network has been set yet.
	 */
	public CyNetwork getReferenceNetwork(){
		if (this.referenceRow < 0){
			return null;
		} else {
			return networkEntries.get(referenceRow).getNetwork();
		}
	}
}
