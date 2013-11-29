package be.svlandeg.diffany.cytoscape.gui;

import java.util.ArrayList;
import java.util.List;

import javax.swing.AbstractListModel;
import javax.swing.ComboBoxModel;
import javax.swing.JComboBox;

import org.cytoscape.model.subnetwork.CyRootNetwork;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.NetworkEntry;

/**
 * Model for the {@link JComboBox} that shows the list of available network collections (or {@link CyRootNetwork}) in 
 * this session. 
 * 
 * @author thpar
 *
 */
public class CollectionDropDownModel extends AbstractListModel implements ComboBoxModel{

	private static final long serialVersionUID = 1L;

	private Model model;
	
	private List<NetworkEntry> collectionEntries = new ArrayList<NetworkEntry>();
	
	private NetworkEntry selectedEntry;
	
	/**
	 * Create a new {@link ComboBoxModel} based on the general {@link Model} of this app and refreshes the 
	 * list of network collections (which on creation will probably be empty).
	 * 
	 * @param model
	 */
	public CollectionDropDownModel(Model model) {
		this.model = model;
		refresh();
	}

	/**
	 * Gathers all available network collections from the {@link Model} and wraps them into
	 * {@link NetworkEntry}s to be presented to the {@link JComboBox}.
	 * 
	 * The refresh will make the combo box GUI redraw.
	 */
	public void refresh(){
		int oldSize = collectionEntries.size();
		collectionEntries = new ArrayList<NetworkEntry>();
		for (CyRootNetwork collection : model.getNetworkCollections()){
			NetworkEntry collectionEntry = new NetworkEntry(collection);
			collectionEntries.add(collectionEntry);
		}
		//select the first entry if available
		if (collectionEntries.size() > 0){
			this.selectedEntry = collectionEntries.get(0);
		}
		//let the gui know all entries might have changed
		this.fireContentsChanged(this, 0, oldSize);
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
