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
	

	/**
	 * Column headers
	 */
	private String[] columns = {"Include", "Network", "Reference"};


	private Model model;
	
	/**
	 * Create a new model based on the general {@link Model} and add it as an {@link Observer} the contained
	 * {@link GUIModel}.
	 * @param model
	 */
	public SelectionTableModel(Model model) {
		this.model = model;
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
			if (!entry.isReference()){
				model.getCurrentProject().setReferenceNetwork(entry.getNetwork());
//				this.fireTableDataChanged();
			} 
			break;
		}
	}


	
	/**
	 * Reload the {@link NetworkEntry}s based on the selected network collection.
	 */
	public void refresh() {
		List<CySubNetwork> subNets = model.getSelectedCollection().getSubNetworkList();
		System.out.println("Number of subnets: "+subNets.size());
		networkEntries = new ArrayList<NetworkEntry>();
		
		//for now, we simply create a whole new CyProject
		//should become more flexible
		CyProject newProject = model.resetProject();
		System.out.println("Refreshing network entries 1");
		
		for (CySubNetwork subNet : subNets){
			System.out.println("Subnet adding 1");
			NetworkEntry entry = new NetworkEntry(subNet);
			entry.setSelected(true);
			entry.setReference(false);
			networkEntries.add(entry);
			System.out.println("Subnet adding 2");
			System.out.println("New project"+newProject);
			newProject.addConditionalNetwork(subNet);
			System.out.println("Subnet adding 3");
		}
		System.out.println("Refreshing network entries 2");
		//select first entry as reference
		if (networkEntries.size() > 0){
			networkEntries.get(0).setReference(true);
			newProject.setReferenceNetwork(networkEntries.get(0).getNetwork());
			newProject.removeConditionalNetwork(networkEntries.get(0).getNetwork());
		}
		System.out.println("Refreshing network entries 3");
		this.fireTableDataChanged();
	}
	
	
}
