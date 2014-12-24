package be.svlandeg.diffany.cytoscape.gui;

/*
 * #%L
 * Diffany
 * %%
 * Copyright (C) 2014 PSB/UGent - Sofie Van Landeghem and Thomas Van Parys
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.AbstractListModel;
import javax.swing.ComboBoxModel;
import javax.swing.JComboBox;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;

/**
 * Model for the {@link JComboBox} that shows the list of available {@link CyProject}s in 
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
	 * 
	 * Create a new {@link ComboBoxModel} based on the general {@link Model} of this app and populates (refreshes) the 
	 * list of network collections.
	 * 
	 * @param model the general model
	 */
	public ProjectDropDownModel(Model model) {
		this.model = model;
		this.refresh();
	}

	/**
	 * Refreshes this model, based on the {@link CyProject}s in the {@link Model}.
	 * Also makes sure the selected project exists and notifies the GUI of the changes.
	 */
	public void refresh(){
		int oldSize = this.entries.size();
		//add new projects
		for (CyProject project : model.getProjects()){
			if (!entries.contains(project)){
				entries.add(project);				
			}
		}
		
		//remove deleted projects
		Set<CyProject> toRemove = new HashSet<CyProject>();
		for (CyProject entry : entries){
			if (!model.getProjects().contains(entry)){
				toRemove.add(entry);
			}
		}
		entries.removeAll(toRemove);
		
		this.selectedEntry = model.getSelectedProject();
		
		if (!entries.contains(this.selectedEntry)){
			this.selectedEntry = null;
		}
		
		if (entries.size() == 0){
			empty = true;
		} else {
			empty = false;
			//set an entry as selected
			if (this.selectedEntry == null){
				this.selectedEntry = entries.get(0);
			}
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
	 * @return whether or not there are any collections in the liest
	 */
	public boolean hasEntries(){
		return !empty;
	}

	


}
