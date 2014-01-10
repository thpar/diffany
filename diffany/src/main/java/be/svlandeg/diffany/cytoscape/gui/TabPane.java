package be.svlandeg.diffany.cytoscape.gui;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Observable;
import java.util.Observer;

import javax.swing.BoxLayout;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;

import org.cytoscape.application.swing.CytoPanelComponent;
import org.cytoscape.application.swing.CytoPanelName;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.swing.DialogTaskManager;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.CyProject.ComparisonMode;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.NetworkEntry;
import be.svlandeg.diffany.cytoscape.actions.RunProjectAction;
import be.svlandeg.diffany.cytoscape.tasks.RunProjectTaskFactory;
import be.svlandeg.diffany.cytoscape.tasks.UpdateVisualStyleTaskFactory;

/**
 * The Control Panel for Diffany (left Cytoscape tab).
 *  
 * @author thpar
 *
 */
public class TabPane extends JPanel implements CytoPanelComponent, Observer, ActionListener, TableModelListener, ChangeListener{

	private static final long serialVersionUID = 1L;
	private Model model;
	private JComboBox collectionDropDown;
	private CollectionDropDownModel comboModel;
	private SelectionTableModel selectionModel;

	private final String COLLECTION_ACTION = "collection";
	private final String MODE_ACTION = "mode";
	private final String RUN_ACTION = "run";
	private JButton runButton;
	
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
		this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		this.add(createCollectionSelectionPanel());
		this.add(createNetworkSelectionPanel());
		this.add(createOptionPanel());
		
		runButton = new JButton(new RunProjectAction(model));
//		runButton.setActionCommand(RUN_ACTION);
//		runButton.addActionListener(this);
		this.add(runButton);
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
		
		collectionDropDown.setActionCommand(COLLECTION_ACTION);
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
	
	/**
	 * Creates the panel containing extra project options like cut off and comparison mode
	 * @return the panel to be added to the side pane
	 */
	private Component createOptionPanel() {
		JPanel panel = new JPanel();
		
		JPanel modePanel = new JPanel();
		JLabel label = new JLabel("Comparison mode: ");
		modePanel.add(label);		
		
		ModeDropDownModel modeDropDownModel = new ModeDropDownModel();
		JComboBox modeDropDown = new JComboBox(modeDropDownModel);
		
		modeDropDown.setSelectedItem(model.getCurrentProject().getMode());
		modeDropDown.setActionCommand(MODE_ACTION);
		modeDropDown.addActionListener(this);
		modePanel.add(modeDropDown);
		
		panel.add(modePanel);
		
		JPanel cutoffPanel = new JPanel();
		SpinnerNumberModel spinModel = new SpinnerNumberModel(model.getCurrentProject().getCutoff(),
				0, 1000, 0.01);
		JSpinner cutoffSpinner = new JSpinner(spinModel);
		cutoffPanel.add(new JLabel("Cutoff"));
		cutoffPanel.add(cutoffSpinner);
		cutoffSpinner.addChangeListener(this);
		
		panel.add(cutoffPanel);
		
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
		this.refreshCyProject();
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		String action = e.getActionCommand();
		if (action.equals(COLLECTION_ACTION)){
			//triggered on collection dropdown action
			JComboBox source = (JComboBox)e.getSource();
			NetworkEntry entry = (NetworkEntry)source.getSelectedItem();
			model.setSelectedCollection((CyRootNetwork)entry.getNetwork());			
		} else if (action.equals(RUN_ACTION)){
			//triggered on Start button click
			
//			RunProjectTaskFactory tf = new RunProjectTaskFactory(model);
//			
//			if (tf.isReady()){
//				TaskIterator it = tf.createTaskIterator();			
//				DialogTaskManager dtm = model.getServices().getDialogTaskManager();
//				dtm.execute(it);
//			}
		} else if (action.equals(MODE_ACTION)){
			JComboBox source = (JComboBox)e.getSource();
			ComparisonMode mode = (CyProject.ComparisonMode)source.getSelectedItem();
			model.getCurrentProject().setMode(mode);
		}
	}

	@Override
	public void tableChanged(TableModelEvent e) {
		// triggered when the data of the table has changed
		this.refreshCyProject();
	}
	
	public void refreshCyProject(){
		CyProject cyProject = model.getCurrentProject();
		cyProject.setConditionalNetworks(this.selectionModel.getConditionalNetworks());
		cyProject.setReferenceNetwork(this.selectionModel.getReferenceNetwork());
		this.runButton.setEnabled(cyProject.canExecute());
		
		UpdateVisualStyleTaskFactory tf = new UpdateVisualStyleTaskFactory(model, cyProject);
		TaskIterator it = tf.createTaskIterator();
		DialogTaskManager dtm = model.getServices().getDialogTaskManager();
		dtm.execute(it);		
	}

	@Override
	public void stateChanged(ChangeEvent e) {
		//triggered when cutoff spinner is changed
		JSpinner spinner = (JSpinner)e.getSource();
		this.model.getCurrentProject().setCutoff((Double)spinner.getValue());
		
	}

	
	
}
