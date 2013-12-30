package be.svlandeg.diffany.cytoscape.gui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Observable;
import java.util.Observer;

import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;

import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.application.swing.CytoPanelName;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.swing.DialogTaskManager;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.NetworkEntry;
import be.svlandeg.diffany.cytoscape.tasks.RunProjectTaskFactory;

/**
 * The Control Panel for Diffany (left Cytoscape tab).
 *  
 * @author thpar
 *
 */
public class TabPane extends JPanel implements CytoPanelComponent, Observer, ActionListener, TableModelListener{

	private static final long serialVersionUID = 1L;
	private Model model;
	private JComboBox collectionDropDown;
	private CollectionDropDownModel comboModel;
	private SelectionTableModel selectionModel;

	/**
	 * Create {@link JPanel} and register as {@link Observer} for the models.
	 * @param model
	 */
	public TabPane(Model model){
		this.model = model;
		model.addObserver(this);			
		
		createTabPaneContent();
	}
	
	private void createTabPaneContent(){
		this.setLayout(new BorderLayout());
		this.add(createCollectionSelectionPanel(), BorderLayout.NORTH);
		this.add(createNetworkSelectionPanel(), BorderLayout.CENTER);
		
		JButton runButton = new JButton("Start");
		runButton.setActionCommand("run");
		runButton.addActionListener(this);
		this.add(runButton, BorderLayout.SOUTH);
	}
	
	/**
	 * Creates the panel containing the drop down list to select a network collection. 
	 * @return {@link JPanel} containing the dropdown list.
	 */
	private Component createCollectionSelectionPanel() {
		comboModel = new CollectionDropDownModel();
		comboModel.refresh(model.getNetworkCollections());
		collectionDropDown = new JComboBox(comboModel);
		collectionDropDown.setEnabled(comboModel.hasEntries());
		
		JPanel collPanel = new JPanel();
		collPanel.add(new JLabel("Network collection: "));
		collPanel.add(collectionDropDown);
		
		collectionDropDown.setActionCommand("collection");
		collectionDropDown.addActionListener(this);
		return collPanel;
	}
	
	/**
	 * Creates the panel containing the network list
	 * @return {@link JPanel} containing the network list.
	 */
	private Component createNetworkSelectionPanel(){
		JPanel panel = new JPanel();
		selectionModel = new SelectionTableModel();
		if (model.getSelectedCollection() !=null){
			selectionModel.refresh(model.getSelectedCollection().getSubNetworkList());			
		}
		JTable table = new JTable(selectionModel);
		table.setPreferredScrollableViewportSize(new Dimension(300, 400));
		
		table.getColumnModel().getColumn(0).setPreferredWidth(20);
		table.getColumnModel().getColumn(2).setPreferredWidth(20);
		table.setFillsViewportHeight(true);
		
		selectionModel.addTableModelListener(this);
		
		JScrollPane scrollPane = new JScrollPane(table);
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
		//triggered on model update
		this.comboModel.refresh(model.getNetworkCollections());
		this.selectionModel.refresh(model.getSelectedCollection().getSubNetworkList());
		this.collectionDropDown.setEnabled(comboModel.hasEntries());
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		String action = e.getActionCommand();
		if (action.equals("collection")){
			//triggered on collection dropdown action
			System.out.println("ComboBox dropdown action");
			JComboBox source = (JComboBox)e.getSource();
			NetworkEntry entry = (NetworkEntry)source.getSelectedItem();
			model.setSelectedCollection((CyRootNetwork)entry.getNetwork());			
		} else if (action.equals("run")){
			//triggered on Start button click
			RunProjectTaskFactory tf = new RunProjectTaskFactory(model);
			
			if (tf.isReady()){
				TaskIterator it = tf.createTaskIterator();			
				DialogTaskManager dtm = model.getServices().getDialogTaskManager();
				dtm.execute(it);
			}
		}
	}

	@Override
	public void tableChanged(TableModelEvent e) {
		// triggered when the data of the table has changed
		CyProject cyProject = model.getCurrentProject();
		cyProject.setConditionalNetworks(this.selectionModel.getConditionalNetworks());
		cyProject.setReferenceNetwork(this.selectionModel.getReferenceNetwork());
	}

	
	
}
