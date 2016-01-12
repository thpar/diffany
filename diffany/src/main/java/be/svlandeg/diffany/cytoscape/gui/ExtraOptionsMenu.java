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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Observable;
import java.util.Observer;

import javax.swing.ButtonGroup;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JRadioButtonMenuItem;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.Model.OverlapOperator;
import be.svlandeg.diffany.cytoscape.actions.UpdateVisualStyleAction;

/**
 * Menu to be placed under Apps, Diffany, Extra options
 * 
 * @author Thomas Van Parys
 *
 */
public class ExtraOptionsMenu extends JMenu implements ActionListener, Observer{


	private static final long serialVersionUID = 1L;
	private ButtonGroup overlapOperatorGroup;
	private JRadioButtonMenuItem minOperator;
	private JRadioButtonMenuItem maxOperator;
	
	/**
	 * The model of the Diffany App
	 */
	private Model appModel;
	private JMenuItem updateVizItem;
	
	/**
	 * Create menu and menu items
	 * 
	 * @param model the model
	 */
	public ExtraOptionsMenu(Model model){
		super("Extra options");
		this.appModel = model;
		this.appModel.addObserver(this);
		
		//Overlap operator
		JMenu overlapMenu = new JMenu("Score operator");
		overlapOperatorGroup = new ButtonGroup();
		minOperator = new JRadioButtonMenuItem("MIN");
		minOperator.setActionCommand("MIN");
		maxOperator = new JRadioButtonMenuItem("MAX");
		maxOperator.setActionCommand("MAX");
		overlapOperatorGroup.add(minOperator);
		overlapOperatorGroup.add(maxOperator);
		
		minOperator.addActionListener(this);
		maxOperator.addActionListener(this);
		
		
		updateOverlapOperatorMenu();
		
		overlapMenu.add(minOperator);
		overlapMenu.add(maxOperator);
		this.add(overlapMenu);
		
		this.addSeparator();
		
		updateVizItem = new JMenuItem(new UpdateVisualStyleAction(model));
		this.add(updateVizItem);
		updateUpdateVizMenuItem();
		
	}
	
	/**
	 * Update states of overlap option checkboxes depending on model
	 */
	private void updateOverlapOperatorMenu(){
		switch(appModel.getOverlapOperator()){
		case MIN:
			overlapOperatorGroup.setSelected(minOperator.getModel(), true);
			break;
		case MAX:
			overlapOperatorGroup.setSelected(maxOperator.getModel(), true);
			break;
		}
	}
	
	/**
	 * Update state of Update Vizualisation menu item
	 */
	private void updateUpdateVizMenuItem(){
		CyProject selectedProject = this.appModel.getSelectedProject();
		this.updateVizItem.setEnabled(selectedProject != null);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		//for now, triggered on Overlap operator action
		JRadioButtonMenuItem radioButton = (JRadioButtonMenuItem)e.getSource();
		if (radioButton.getActionCommand().equals("MIN")){
			appModel.setOverlapOperator(OverlapOperator.MIN);
		} else {
			appModel.setOverlapOperator(OverlapOperator.MAX);
		}
		
	}

	@Override
	public void update(Observable o, Object arg) {
		//triggers on model update
		this.updateUpdateVizMenuItem();
	}

}
