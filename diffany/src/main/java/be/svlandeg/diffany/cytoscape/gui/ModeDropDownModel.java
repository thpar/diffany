package be.svlandeg.diffany.cytoscape.gui;

import javax.swing.AbstractListModel;
import javax.swing.ComboBoxModel;
import javax.swing.JComboBox;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;

/**
 * Model for the {@link JComboBox} displaying available algorithm modes.
 *  
 * @author Thomas Van Parys
 *
 */
public class ModeDropDownModel extends AbstractListModel implements ComboBoxModel {

	
	private static final long serialVersionUID = 1L;
	private Model.ComparisonMode selected;
	
	@Override
	public int getSize() {
		return Model.ComparisonMode.values().length;
	}

	@Override
	public Object getElementAt(int index) {
		return Model.ComparisonMode.values()[index];
	}

	@Override
	public void setSelectedItem(Object anItem) {
		this.selected = (Model.ComparisonMode)anItem;
	}

	@Override
	public Object getSelectedItem() {
		return selected;
	}

}
