package be.svlandeg.diffany.cytoscape.gui;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import javax.swing.AbstractListModel;
import javax.swing.ComboBoxModel;
import javax.swing.JComboBox;

import org.cytoscape.model.subnetwork.CyRootNetwork;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;

/**
 * Model for the {@link JComboBox} that shows the list of available network collections (or {@link CyRootNetwork}) in 
 * this session. 
 * 
 * @author Thomas Van Parys
 *
 */
public class ProjectDropDownModel extends AbstractListModel implements ComboBoxModel{

	private static final long serialVersionUID = 1L;

	private static final String EMPTY_MESSAGE = "No Collections";
	
	private boolean empty = true;
	
	private List<CyProject> entries = new ArrayList<CyProject>();
	
	private CyProject selectedEntry;

	private Model model;
	
	/**
	 * TODO: update documentation
	 * 
	 * Create a new {@link ComboBoxModel} based on the general {@link Model} of this app and refreshes the 
	 * list of network collections (which on creation will probably be empty).
	 * 
	 */
	public ProjectDropDownModel(Model model) {
		this.model = model;
	}

	public void refresh(){
		int oldSize = this.entries.size();
		for (CyProject project : model.getProjects()){
			entries.add(project);
		}
		
		if (model.getProjects().size() == 0){
			empty = true;
		} else {
			empty = false;
		}
		this.fireContentsChanged(this, 0, oldSize);
	}
	
	@Override
	public int getSize() {
		if (empty){
			//make a spot for the "Empty" entry
			return 1;
		} else {
			return entries.size();
		}
	}

	@Override
	public Object getElementAt(int index) {
		if (empty){
			return EMPTY_MESSAGE;
		} else {
			return entries.get(index);			
		}
	}

	@Override
	public void setSelectedItem(Object anItem) {
		this.selectedEntry = (CyProject)anItem;		
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
