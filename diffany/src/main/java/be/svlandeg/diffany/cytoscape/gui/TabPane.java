package be.svlandeg.diffany.cytoscape.gui;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Observable;
import java.util.Observer;
import java.util.Set;

import javax.swing.Icon;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;

import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.application.swing.CytoPanelName;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;

import be.svlandeg.diffany.NetworkEntry;
import be.svlandeg.diffany.cytoscape.Model;

public class TabPane extends JPanel implements CytoPanelComponent, Observer, ActionListener{

	private static final long serialVersionUID = 1L;
	private Model model;
	private GUIModel guiModel;

	public TabPane(Model model, GUIModel guiModel){
		this.model = model;
		this.guiModel = guiModel;
		model.addObserver(this);
		guiModel.addObserver(this);
			
		
		this.add(createCollectionSelectionPanel());

		this.add(createNetworkSelectionPanel());
	}
	
	private Component createCollectionSelectionPanel() {
		JPanel panel = new JPanel();
		JComboBox dropDown = new JComboBox();
		repopulate(dropDown);
		panel.add(dropDown);
		
		dropDown.addActionListener(this);
		
		return panel;
	}

	private Component createNetworkSelectionPanel(){
		JPanel panel = new JPanel();
		JTable table = new JTable(new SelectionTableModel(guiModel));
		JScrollPane scrollPane = new JScrollPane(table);
		table.setFillsViewportHeight(true);
		panel.add(scrollPane);
		return panel;
	}
	
	
	
	private void repopulate(JComboBox dropDown) {
		Set<CyRootNetwork> rootNets = model.getNetworkCollections();
		for (CyRootNetwork root : rootNets){
			NetworkEntry entry = new NetworkEntry(root);
			dropDown.addItem(entry);
		}
	}

	@Override
	public Component getComponent() {
		return this;
	}

	@Override
	public CytoPanelName getCytoPanelName() {
		return CytoPanelName.WEST;
	}

	@Override
	public String getTitle() {
		return new String("Diffany");
	}

	@Override
	public Icon getIcon() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void update(Observable o, Object arg) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		JComboBox source = (JComboBox)e.getSource();
		NetworkEntry entry = (NetworkEntry)source.getSelectedItem();
		guiModel.setSelectedCollection(entry.getNetwork());
	}
	
	
}
