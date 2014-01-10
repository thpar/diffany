package be.svlandeg.diffany.cytoscape.gui;

import javax.swing.AbstractListModel;
import javax.swing.ComboBoxModel;
import javax.swing.JComboBox;

import be.svlandeg.diffany.cytoscape.CyProject;

/**
 * Model for the {@link JComboBox} displaying available alorithm modes.
 *  
 * @author Thomas Van Parys
 *
 */
public class ModeDropDownModel extends AbstractListModel implements ComboBoxModel {

	
	private static final long serialVersionUID = 1L;
	private CyProject.ComparisonMode selected;
	
	@Override
	public int getSize() {
		return CyProject.ComparisonMode.values().length;
	}

	@Override
	public Object getElementAt(int index) {
		return CyProject.ComparisonMode.values()[index];
	}

	@Override
	public void setSelectedItem(Object anItem) {
		this.selected = (CyProject.ComparisonMode)anItem;
	}

	@Override
	public Object getSelectedItem() {
		return selected;
	}

}
