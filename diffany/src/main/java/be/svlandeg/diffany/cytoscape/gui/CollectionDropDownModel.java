package be.svlandeg.diffany.cytoscape.gui;

import java.util.ArrayList;
import java.util.List;

import javax.swing.AbstractListModel;
import javax.swing.ComboBoxModel;

import org.cytoscape.model.subnetwork.CyRootNetwork;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.NetworkEntry;

public class CollectionDropDownModel extends AbstractListModel implements ComboBoxModel{


	private static final long serialVersionUID = 1L;

	private Model model;
	
	private List<NetworkEntry> collectionEntries = new ArrayList<NetworkEntry>();
	
	private NetworkEntry selectedEntry;
	
	public CollectionDropDownModel(Model model) {
		this.model = model;
		refresh();
	}

	public void refresh(){
		System.out.println("refresh in");
		int oldSize = collectionEntries.size();
		collectionEntries = new ArrayList<NetworkEntry>();
		for (CyRootNetwork collection : model.getNetworkCollections()){
			NetworkEntry collectionEntry = new NetworkEntry(collection);
			collectionEntries.add(collectionEntry);
		}
		System.out.println("fire in the hole!");
		//FIXME this is going wrong
		this.fireContentsChanged(this, 0, oldSize);
		System.out.println("refresh out");
	}
	
	@Override
	public int getSize() {
		int size = collectionEntries.size();
		return size;
	}

	@Override
	public Object getElementAt(int index) {
		return collectionEntries.get(index);
	}

	@Override
	public void setSelectedItem(Object anItem) {
		this.selectedEntry = (NetworkEntry)anItem;		
	}

	@Override
	public Object getSelectedItem() {
		return this.selectedEntry;
	}


}
