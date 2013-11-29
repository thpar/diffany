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

/**
 * The Control Panel for Diffany (left Cytoscape tab).
 *  
 * @author thpar
 *
 */
public class TabPane extends JPanel implements CytoPanelComponent, Observer, ActionListener, NetworkAddedListener{

	private static final long serialVersionUID = 1L;
	private Model model;
	private JComboBox collectionDropDown;
	private CollectionDropDownModel comboModel;

	/**
	 * Create {@link JPanel} and register as {@link Observer} for the models.
	 * @param model
	 */
	public TabPane(Model model){
		this.model = model;
		model.addObserver(this);			
		model.getGuiModel().addObserver(this);
		
		this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		
		this.add(createCollectionSelectionPanel());
		this.add(createNetworkSelectionPanel());
	}
	
	/**
	 * Creates the panel containing the drop down list to select a network collection. 
	 * @return {@link JPanel} containing the dropdown list.
	 */
	private Component createCollectionSelectionPanel() {
		JPanel panel = new JPanel();
		comboModel = new CollectionDropDownModel(model);
		collectionDropDown = new JComboBox(comboModel);
		panel.add(collectionDropDown);
		
		collectionDropDown.addActionListener(this);
		return panel;
	}
	
	/**
	 * Creates the panel containing the network list
	 * @return {@link JPanel} containing the network list.
	 */
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
		//triggered on model or guiModel change
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		//triggered on collection dropdown action
		JComboBox source = (JComboBox)e.getSource();
		NetworkEntry entry = (NetworkEntry)source.getSelectedItem();
		model.getGuiModel().setSelectedCollection(entry.getNetwork());
	}

	@Override
	public void handleEvent(NetworkAddedEvent e) {
		//triggered on Cytoscape NetworkAdded
		comboModel.refresh();
	}
	
}
