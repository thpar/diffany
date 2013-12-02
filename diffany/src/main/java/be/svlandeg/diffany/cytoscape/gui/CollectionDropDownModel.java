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

	private static final Object EMPTY_MESSAGE = "No Collections";

	private Model model;
	
	private List<NetworkEntry> collectionEntries = new ArrayList<NetworkEntry>();
	
	private NetworkEntry selectedEntry;
	
	private boolean empty = true;
	
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
		int oldSize = this.getSize();
		collectionEntries = new ArrayList<NetworkEntry>();
		for (CyRootNetwork collection : model.getNetworkCollections()){
			NetworkEntry collectionEntry = new NetworkEntry(collection);
			collectionEntries.add(collectionEntry);
		}
		//select the first entry if available
		if (collectionEntries.size() > 0){
			this.empty = false;
			this.selectedEntry = collectionEntries.get(0);
		} else {
			empty = true;
		}
		//let the gui know all entries might have changed
		this.fireContentsChanged(this, 0, oldSize);
	}
	
	@Override
	public int getSize() {
		if (empty){
			return 1;
		} else {
			return collectionEntries.size();
		}
	}

	@Override
	public Object getElementAt(int index) {
		if (empty){
			return EMPTY_MESSAGE;
		} else {
			return collectionEntries.get(index);			
		}
	}

	@Override
	public void setSelectedItem(Object anItem) {
		this.selectedEntry = (NetworkEntry)anItem;		
	}

	@Override
	public Object getSelectedItem() {
		if (empty){
			return EMPTY_MESSAGE;
		}
		return this.selectedEntry;
	}
	
	/**
	 * Checks if there are any collections in the list
	 * @return
	 */
	public boolean hasEntries(){
		return !empty;
	}


}
