package be.svlandeg.diffany.cytoscape.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.ButtonGroup;
import javax.swing.JMenu;
import javax.swing.JRadioButtonMenuItem;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.Model.OverlapOperator;

public class ExtraOptionsMenu extends JMenu implements ActionListener {


	private static final long serialVersionUID = 1L;
	private ButtonGroup overlapOperatorGroup;
	private JRadioButtonMenuItem minOperator;
	private JRadioButtonMenuItem maxOperator;
	private Model appModel;
	
	public ExtraOptionsMenu(Model model){
		super("Extra options");
		this.appModel = model;
		
		//Overlap operator
		JMenu overlapMenu = new JMenu("Overlap operator");
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
		
	}
	
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

}
