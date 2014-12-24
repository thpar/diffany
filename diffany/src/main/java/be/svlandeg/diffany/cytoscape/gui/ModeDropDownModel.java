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

import javax.swing.AbstractListModel;
import javax.swing.ComboBoxModel;
import javax.swing.JComboBox;

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
