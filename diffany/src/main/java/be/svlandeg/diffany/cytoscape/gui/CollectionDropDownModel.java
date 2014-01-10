package be.svlandeg.diffany.cytoscape.gui;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import javax.swing.AbstractListModel;
import javax.swing.ComboBoxModel;
import javax.swing.JComboBox;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.NetworkEntry;

/**
 * Model for the {@link JComboBox} that shows the list of available network collections (or {@link CyRootNetwork}) in 
 * this session. 
 * 
 * @author Thomas Van Parys
 *
 */
public class CollectionDropDownModel extends AbstractListModel implements ComboBoxModel{

	private static final long serialVersionUID = 1L;

	private static final Object EMPTY_MESSAGE = "No Collections";

	private List<NetworkEntry> collectionEntries = new ArrayList<NetworkEntry>();
	
	private NetworkEntry selectedEntry;
	
	private boolean empty = true;
	
	/**
	 * Create a new {@link ComboBoxModel} based on the general {@link Model} of this app and refreshes the 
	 * list of network collections (which on creation will probably be empty).
	 * 
	 */
	public CollectionDropDownModel() {
	}

	/**
	 * Takes a set of {@link CyRootNetwork}s and refreshes the entries in this ComboBox.
	 * A previously selected entry will stay selected if the corresponding network still exists.
	 * 
	 * The refresh will make the combo box GUI redraw.
	 * 
	 * @param collections the network collections currently available in this Cytoscape session.
	 */
	public void refresh(Set<CyRootNetwork> collections){
		int oldSize = this.getSize();
		List<NetworkEntry> newCollectionEntries = new ArrayList<NetworkEntry>();
		for (CyRootNetwork collection : collections){
			NetworkEntry entry = getEntry(collection);
			if (entry != null){
				newCollectionEntries.add(entry);
			} else {
				NetworkEntry newEntry = new NetworkEntry(collection);
				newCollectionEntries.add(newEntry);				
			}
		}
		this.collectionEntries = newCollectionEntries;
		//select the first entry if available, only if the previously selected entry was not available
		if (selectedEntry == null || !this.collectionEntries.contains(this.selectedEntry)){
			if (collectionEntries.size() > 0){
				this.empty = false;
				this.selectedEntry = collectionEntries.get(0);
			}
		}
		this.empty = this.collectionEntries.size()==0;
		
		//let the gui know all entries might have changed
		this.fireContentsChanged(this, 0, oldSize);
	}
	
	private NetworkEntry getEntry(CyNetwork network){
		for (NetworkEntry entry : this.collectionEntries){
			if (network == entry.getNetwork()){
				return entry;
			}
		}
		return null;
	}
	
	@Override
	public int getSize() {
		if (empty){
			//make a spot for the "Empty" entry
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
