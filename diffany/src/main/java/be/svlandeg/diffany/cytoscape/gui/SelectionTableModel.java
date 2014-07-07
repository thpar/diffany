package be.svlandeg.diffany.cytoscape.gui;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.table.AbstractTableModel;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CySubNetwork;

import be.svlandeg.diffany.cytoscape.CyProject;

/**
 * The model for the table of available networks in the selected collection.
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
	
	private int referenceRow = -1;

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
		switch(columnIndex){
		case 1:
			return false;
		case 0:
		case 2:
			return true;
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
			//a reference network should be selected anyway
			if (!entry.isSelected() && check){
				entry.setSelected(true);
			}
			
			if (check){
				this.setReference(rowIndex);				
			} else {
				this.unselectReference(rowIndex);
			}
			break;
		}
	}



	/**
	 * Set the given row to contain the reference network
	 * @param row row number of reference network
	 */
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
	
	private void unselectReference(int row) {
		this.referenceRow = -1;
		this.networkEntries.get(row).setReference(false);
		this.fireTableDataChanged();
	}

	/**
	 * Reload the {@link NetworkEntry}s based on the subnetworks from the selected network collections.
	 * Compare with the previous list to maintain selections.
	 * 
	 * @param list of {@link CySubNetwork}s to be displayed.
	 */
	public void refresh(CyProject project) {
		List<CySubNetwork> subNets = project.getCollection().getSubNetworkList();
		
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
					entry.setSelected(
							project.isConditionalNetwork(subNet) ||
							project.isReferenceNetwork(subNet));
					networkEntries.add(entry);
					entryMap.put(subNet, entry);
					if (project.isReferenceNetwork(subNet)){
						this.referenceRow = networkEntries.indexOf(entry); 
						entry.setReference(true);
					} else {
						entry.setReference(false);
					}
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
	 * @return the {@link CyNetwork} to be used as reference network. Returns null if no reference network is set.
	 */
	public CyNetwork getReferenceNetwork(){
		if (this.referenceRow < 0){
			return null;
		} else {
			return networkEntries.get(referenceRow).getNetwork();
		}
	}

	/**
	 * Clears all data in the model. To be used when no network collection could be selected.
	 */
	public void clear() {
		this.networkEntries = new ArrayList<NetworkEntry>();
		this.entryMap = new HashMap<CyNetwork, NetworkEntry>();
		this.referenceRow = -1;
		this.fireTableDataChanged();
	}
	
	/**
	 * Get the entry at given row
	 * 
	 * @param row
	 * @return 
	 */
	public NetworkEntry getNetworkEntry(int row){
		return this.networkEntries.get(row);		
	}
	
	/**
	 * Get the row number of a given {@link CyNetwork}
	 * 
	 * @param net the {@link CyNetwork} to retrieve the row number of
	 * @return row number of given {@link CyNetwork}. Returns -1 if no network was found.
	 */
	public int getRowNumber(CyNetwork net){
		NetworkEntry entry = entryMap.get(net);
		if (entry !=null){
			return networkEntries.indexOf(entry);			
		} else {
			return -1;
		}
	}
}
