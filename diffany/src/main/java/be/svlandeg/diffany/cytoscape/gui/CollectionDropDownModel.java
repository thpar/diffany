package be.svlandeg.diffany.cytoscape.gui;

import java.util.ArrayList;
import java.util.List;

import javax.swing.DefaultComboBoxModel;
import javax.swing.event.ListDataListener;

import org.cytoscape.model.subnetwork.CyRootNetwork;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.NetworkEntry;

public class CollectionDropDownModel extends DefaultComboBoxModel {


	private static final long serialVersionUID = 1L;

	private Model model;
	
	private List<NetworkEntry> collectionEntries;
	
	public CollectionDropDownModel(Model model) {
		this.model = model;
		refresh();
	}

	public void refresh(){
		collectionEntries = new ArrayList<NetworkEntry>();
		for (CyRootNetwork collection : model.getNetworkCollections()){
			NetworkEntry collectionEntry = new NetworkEntry(collection);
			collectionEntries.add(collectionEntry);
		}
	}
	
	@Override
	public int getSize() {
		return collectionEntries.size();
	}

	@Override
	public Object getElementAt(int index) {
		return collectionEntries.get(index);
	}



}
