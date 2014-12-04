package be.svlandeg.diffany.cytoscape.gui;

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
 * Menu to be placed under Apps > Diffany > Extra options
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
	 * @param model
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
