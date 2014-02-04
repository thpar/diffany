package be.svlandeg.diffany.cytoscape.gui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Observable;
import java.util.Observer;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
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
import be.svlandeg.diffany.cytoscape.actions.UpdateVisualStyleAction;
import be.svlandeg.diffany.cytoscape.tasks.UpdateVisualStyleTaskFactory;

/**
 * The Control Panel for Diffany (left Cytoscape tab).
 *  
 * @author Thomas Van Parys
 *
 */
public class TabPane extends JPanel implements CytoPanelComponent, Observer, ActionListener, 
									TableModelListener, ChangeListener, ListSelectionListener{

	private static final long serialVersionUID = 1L;
	private Model model;
	private JComboBox collectionDropDown;
	private CollectionDropDownModel collectionModel;
	private SelectionTableModel selectionModel;

	private final String COLLECTION_ACTION = "collection";
	private final String MODE_ACTION = "mode";
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
		
		this.add(createNetworkSelectionPanel());
		
		this.add(createOptionPanel());
		
		JPanel runPanel = new JPanel();
		runPanel.setBorder(BorderFactory.createTitledBorder("Run Diffany"));
		runButton = new JButton(new RunProjectAction(model));
		runButton.setEnabled(model.getCurrentProject().canExecute());
		
		runPanel.add(new JButton(new UpdateVisualStyleAction(model, model.getCurrentProject())));
		runPanel.add(runButton);
		
		
		this.add(runPanel);
	}
	

	
	
	/**
	 * Creates the panel containing the network list
	 * @return {@link JPanel} containing the network list.
	 */
	private Component createNetworkSelectionPanel(){
		JPanel panel = new JPanel();
		panel.setLayout(new BorderLayout());
		panel.setBorder(BorderFactory.createTitledBorder("Input networks"));
		
		collectionModel = new CollectionDropDownModel();
		collectionModel.refresh(model.getNetworkCollections());
		collectionDropDown = new JComboBox(collectionModel);
		collectionDropDown.setEnabled(collectionModel.hasEntries());
		
		JPanel collPanel = new JPanel();
		collPanel.setLayout(new BoxLayout(collPanel, BoxLayout.LINE_AXIS));
		collPanel.add(new JLabel("Network collection: "));
		collPanel.add(collectionDropDown);
		
		collectionDropDown.setActionCommand(COLLECTION_ACTION);
		collectionDropDown.addActionListener(this);
		
		panel.add(collPanel, BorderLayout.NORTH);
		
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
		table.getSelectionModel().addListSelectionListener(this);
		
		JScrollPane scrollPane = new JScrollPane(table);
		panel.add(scrollPane, BorderLayout.CENTER);
		return panel;
	}
	
	/**
	 * Creates the panel containing extra project options like cut off and comparison mode
	 * @return the panel to be added to the side pane
	 */
	private Component createOptionPanel() {
		JPanel panel = new JPanel();
		panel.setLayout(new BoxLayout(panel, BoxLayout.LINE_AXIS));
		panel.setBorder(BorderFactory.createTitledBorder("Options"));
		
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
		this.collectionModel.refresh(model.getNetworkCollections());
		this.collectionDropDown.setEnabled(collectionModel.hasEntries());
		
		if (model.getSelectedCollection() !=null){
			this.selectionModel.refresh(model.getSelectedCollection().getSubNetworkList());			
		} else {
			this.selectionModel.clear();
		}
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		String action = e.getActionCommand();
		if (action.equals(COLLECTION_ACTION)){
			//triggered on collection dropdown action
			JComboBox source = (JComboBox)e.getSource();
			Object selected = source.getSelectedItem();
			if (selected instanceof NetworkEntry){
				NetworkEntry entry = (NetworkEntry)selected;
				model.setSelectedCollection((CyRootNetwork)entry.getNetwork());							
			} else {
				model.setSelectedCollection(null);
			}
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
	
	/**
	 * Reads the selections in the GUI and reflects them in the current {@link CyProject} object.
	 */
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

	@Override
	public void valueChanged(ListSelectionEvent e) {
		ListSelectionModel lsm = (ListSelectionModel)e.getSource();
		
		
	}

	
	
}
