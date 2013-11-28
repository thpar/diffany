package be.svlandeg.diffany.cytoscape.gui;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Observable;
import java.util.Observer;
import java.util.Set;

import javax.swing.BoxLayout;
import javax.swing.ComboBoxModel;
import javax.swing.Icon;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;

import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.application.swing.CytoPanelName;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.events.NetworkAddedEvent;
import org.cytoscape.model.events.NetworkAddedListener;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CyRootNetworkManager;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.NetworkEntry;

public class TabPane extends JPanel implements CytoPanelComponent, Observer, ActionListener, NetworkAddedListener{

	private static final long serialVersionUID = 1L;
	private Model model;
	private JComboBox collectionDropDown;
	private CollectionDropDownModel comboModel;

	public TabPane(Model model){
		this.model = model;
		model.addObserver(this);			
		model.getGuiModel().addObserver(this);
		
		this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		
		this.add(createCollectionSelectionPanel());
		this.add(createNetworkSelectionPanel());
	}
	
	private Component createCollectionSelectionPanel() {
		JPanel panel = new JPanel();
		comboModel = new CollectionDropDownModel(model);
		collectionDropDown = new JComboBox();
		panel.add(collectionDropDown);
		
		collectionDropDown.addActionListener(this);
		return panel;
	}
	

	private Component createNetworkSelectionPanel(){
		JPanel panel = new JPanel();
		JTable table = new JTable(new SelectionTableModel(model));
		JScrollPane scrollPane = new JScrollPane(table);
		table.setFillsViewportHeight(true);
		panel.add(scrollPane);
		return panel;
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
		model.getGuiModel().setSelectedCollection(entry.getNetwork());
	}

	@Override
	public void handleEvent(NetworkAddedEvent e) {
		System.out.println(">>>Network was added!!!");
		comboModel.refresh();
		collectionDropDown.updateUI();
		System.out.println("Network was added!!!<<<");
	}
	
}
