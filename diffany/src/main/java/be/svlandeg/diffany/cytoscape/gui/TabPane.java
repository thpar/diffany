package be.svlandeg.diffany.cytoscape.gui;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collection;
import java.util.HashSet;
import java.util.Observable;
import java.util.Observer;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
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
import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyEdge.Type;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyRow;
import org.cytoscape.model.CyTable;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.view.model.View;
import org.cytoscape.view.presentation.property.BasicVisualLexicon;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.swing.DialogTaskManager;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.Model.ComparisonMode;
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
	private ProjectDropDownModel collectionModel;
	private SelectionTableModel networkTableModel;

	private static final String COLLECTION_ACTION = "collection";
	private static final String MODE_ACTION = "mode";
	private static final String GENERATE_DIFF_ACTION = "generate differential networks";
	private static final String GENERATE_CONSENSUS_ACTION = "generate consensus networks";

	private static final String REQUIRE_REF_NET_ACTION = "require reference network";
	
	private JButton runButton;
	private JButton updateVizButton;
	private JTable networkTable;
	private JCheckBox requireRefNetCheckBox;
	private JSlider supportSlider;
	private JTable filterEdgesTable;
	private EdgeFilterTableModel filterTableModel;
	
	/**
	 * Create {@link JPanel} and register as {@link Observer} for the model.
	 * @param model
	 */
	public TabPane(Model model){
		this.model = model;
		model.addObserver(this);			
		
		createTabPaneContent();
	}
	
	/**
	 * Create the content of the complete panel and adds it to this side bar. 
	 */
	private void createTabPaneContent(){
		this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		
		this.add(createNetworkSelectionPanel());

		this.add(createFilterPanel());
		
		this.add(createOptionPanel());

		JPanel runPanel = new JPanel();
		runPanel.setBorder(BorderFactory.createTitledBorder("Run Diffany"));
		runPanel.setAlignmentX(CENTER_ALIGNMENT);
		runButton = new JButton(new RunProjectAction(model));
		updateVizButton = new JButton(new UpdateVisualStyleAction(model));
		
		runPanel.add(updateVizButton);
		runPanel.add(runButton);
		
		CyProject selectedProject = model.getSelectedProject();
		if (selectedProject!=null){
			runButton.setEnabled(model.getSelectedProject().canExecute(model));
			updateVizButton.setEnabled(true);
		} else {
			runButton.setEnabled(false);
			updateVizButton.setEnabled(false);
		}
		
		this.add(runPanel);
	}
	

	
	
	

	/**
	 * Creates the panel containing the dropdown list with {@link CyProject}s and the table
	 * containing the networks.
	 * 
	 * @return {@link JPanel} containing the project and network lists.
	 */
	private Component createNetworkSelectionPanel(){

		JPanel panel = new JPanel();
		panel.setLayout(new BorderLayout());
		panel.setBorder(BorderFactory.createTitledBorder("Input networks"));
		panel.setAlignmentX(CENTER_ALIGNMENT);
		
		collectionModel = new ProjectDropDownModel(model);
		
		collectionDropDown = new JComboBox(collectionModel);
		collectionDropDown.setEnabled(collectionModel.hasEntries());
		
		JPanel collPanel = new JPanel();
		collPanel.setLayout(new BoxLayout(collPanel, BoxLayout.LINE_AXIS));
		collPanel.add(new JLabel("Network collection: "));
		collPanel.add(collectionDropDown);
		
		collectionDropDown.setActionCommand(COLLECTION_ACTION);
		collectionDropDown.addActionListener(this);
		
		panel.add(collPanel, BorderLayout.NORTH);
		
		networkTableModel = new SelectionTableModel();
		if (model.getSelectedProject() !=null){
			networkTableModel.refresh(model.getSelectedProject());			
		}
		networkTable = new JTable(networkTableModel);
		networkTable.setPreferredScrollableViewportSize(new Dimension(300, 400));
		
		networkTable.getColumnModel().getColumn(0).setPreferredWidth(20);
		networkTable.getColumnModel().getColumn(2).setPreferredWidth(20);
		networkTable.setFillsViewportHeight(true);
		
		networkTableModel.addTableModelListener(this);
		
		networkTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		networkTable.getSelectionModel().addListSelectionListener(this);
		
		CyNetwork focused = model.getNetworkInFocus();
		if (focused !=null){
			int focusRow = networkTableModel.getRowNumber(focused);
			if (focusRow >=0){
				networkTable.setRowSelectionInterval(focusRow, focusRow);				
			}
		}
		
		JScrollPane scrollPane = new JScrollPane(networkTable);
		panel.add(scrollPane, BorderLayout.CENTER);
		return panel;
	}
	
	private Component createFilterPanel() {
		JPanel panel = new JPanel();
		panel.setLayout(new BorderLayout());
		panel.setBorder(BorderFactory.createTitledBorder("Interaction filter"));
		panel.setAlignmentX(CENTER_ALIGNMENT);
		
		filterTableModel = new EdgeFilterTableModel();
		if (model.getSelectedProject() != null){
			filterTableModel.refresh(model.getSelectedProject());
		}
		
		filterEdgesTable = new JTable(filterTableModel);
		filterEdgesTable.setPreferredScrollableViewportSize(new Dimension(300, 100));
		filterEdgesTable.getColumnModel().getColumn(1).setPreferredWidth(20);
		filterEdgesTable.setFillsViewportHeight(true);
		
		filterTableModel.addTableModelListener(this);
		
		filterEdgesTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		//table.getSelectionModel().addListSelectionListener(this);
		JScrollPane scrollPane = new JScrollPane(filterEdgesTable);
		panel.add(scrollPane, BorderLayout.CENTER);
				
		return panel;
	}
	
	
	/**
	 * Creates the panel containing extra project options like cut off and comparison mode
	 * @return the panel to be added to the side pane
	 */
	private Component createOptionPanel() {
		JPanel optionPanel = new JPanel();
		optionPanel.setLayout(new BoxLayout(optionPanel, BoxLayout.PAGE_AXIS));
		optionPanel.setBorder(BorderFactory.createTitledBorder("Options"));
		optionPanel.setMaximumSize(new Dimension(Short.MAX_VALUE, Short.MAX_VALUE));
		optionPanel.setAlignmentX(CENTER_ALIGNMENT);
		
		JPanel topPanel = new JPanel();
		topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.LINE_AXIS));
		topPanel.setAlignmentX(LEFT_ALIGNMENT);
		topPanel.setMaximumSize(new Dimension(Short.MAX_VALUE, Short.MAX_VALUE));
		
		JPanel modePanel = new JPanel();
		JLabel label = new JLabel("Comparison mode: ");
		modePanel.add(label);		
		
		ModeDropDownModel modeDropDownModel = new ModeDropDownModel();
		JComboBox modeDropDown = new JComboBox(modeDropDownModel);
		
		modeDropDown.setSelectedItem(model.getMode());
		modeDropDown.setActionCommand(MODE_ACTION);
		modeDropDown.addActionListener(this);
		modePanel.add(modeDropDown);
		
		topPanel.add(modePanel);
		
		JPanel cutoffPanel = new JPanel();
		SpinnerNumberModel spinModel = new SpinnerNumberModel(model.getCutoff(),
				0, 1000, 0.01);
		JSpinner cutoffSpinner = new JSpinner(spinModel);
		cutoffPanel.add(new JLabel("Cutoff"));
		cutoffPanel.add(cutoffSpinner);
		cutoffSpinner.addChangeListener(this);
		
		topPanel.add(cutoffPanel);
		optionPanel.add(topPanel);
		
		
		JPanel outputSelectionPanel = new JPanel();
		outputSelectionPanel.setAlignmentX(LEFT_ALIGNMENT);
		outputSelectionPanel.setMaximumSize(new Dimension(Short.MAX_VALUE, Short.MAX_VALUE));
		outputSelectionPanel.setLayout(new BoxLayout(outputSelectionPanel, BoxLayout.PAGE_AXIS));
		
		JCheckBox generateDiffNetCheckBox = new JCheckBox("Differential networks");
		generateDiffNetCheckBox.setActionCommand(GENERATE_DIFF_ACTION);
		generateDiffNetCheckBox.addActionListener(this);
		generateDiffNetCheckBox.setSelected(model.isGenerateDiffNets());
		outputSelectionPanel.add(generateDiffNetCheckBox);
		
		JCheckBox generateConsensusNetCheckBox = new JCheckBox("Consensus networks");
		generateConsensusNetCheckBox.setActionCommand(GENERATE_CONSENSUS_ACTION);
		generateConsensusNetCheckBox.addActionListener(this);
		generateConsensusNetCheckBox.setSelected(model.isGenerateConsensusNets());
		outputSelectionPanel.add(generateConsensusNetCheckBox);
		
		optionPanel.add(outputSelectionPanel);
		
		optionPanel.add(Box.createRigidArea(new Dimension(0,10)));
		JPanel supportPanel = new JPanel();
		supportPanel.setLayout(new BoxLayout(supportPanel, BoxLayout.PAGE_AXIS));
		supportPanel.setAlignmentX(LEFT_ALIGNMENT);
		JLabel supportLabel = new JLabel("Consensus edges cutoff:");
		supportSlider = new JSlider(JSlider.HORIZONTAL, 0,4,4);
		supportSlider.addChangeListener(this);
		supportSlider.setMajorTickSpacing(1);
		supportSlider.setMinorTickSpacing(1);
		supportSlider.setPaintTicks(true);
		supportSlider.setPaintLabels(true);
		supportSlider.setEnabled(false);
		
		requireRefNetCheckBox = new JCheckBox("Require reference network edge");
		requireRefNetCheckBox.setActionCommand(REQUIRE_REF_NET_ACTION);
		requireRefNetCheckBox.addActionListener(this);
		requireRefNetCheckBox.setSelected(true);
		requireRefNetCheckBox.setEnabled(false);
		
		supportPanel.add(supportLabel);
		supportPanel.add(supportSlider);
		supportPanel.add(requireRefNetCheckBox);
		optionPanel.add(supportPanel);
				
		return optionPanel;
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
		this.collectionModel.refresh();
		this.collectionDropDown.setEnabled(collectionModel.hasEntries());
		
		if (model.getSelectedProject() !=null){
			this.networkTableModel.refresh(model.getSelectedProject());
			this.filterTableModel.refresh(model.getSelectedProject());
			
			CyNetwork focused = model.getNetworkInFocus();
			if (focused !=null){
				int focusRow = networkTableModel.getRowNumber(focused);
				if (focusRow >=0){
					networkTable.setRowSelectionInterval(focusRow, focusRow);					
				}
			}
			
			refreshCyProject();
		} else {
			this.networkTableModel.clear();
			this.filterTableModel.clear();
		}
		updateSupportSlider();
	}
	
	private void updateSupportSlider(){
		//enable extra options
		this.supportSlider.setEnabled(model.getSelectedProject()!=null);
		if (model.getSelectedProject()!=null 
				&& !model.isGenerateDiffNets() 
				&& model.getSelectedProject().getReferenceNetwork()!=null){
			this.requireRefNetCheckBox.setEnabled(true);
			this.requireRefNetCheckBox.setSelected(true);
		} else {
			this.requireRefNetCheckBox.setEnabled(false);
			this.requireRefNetCheckBox.setSelected(false);
		}
		
		this.requireRefNetCheckBox.setSelected(model.getSelectedProject().getReferenceNetwork() != null);
		int currentTicks = supportSlider.getMaximum();
		int newNumberOfTicks = this.model.getSelectedProject().getNumberOfInputNetworks();
		if (currentTicks != newNumberOfTicks){
			this.supportSlider.setMaximum(newNumberOfTicks);
			this.supportSlider.setMinimum(newNumberOfTicks/2);
			this.supportSlider.setValue(newNumberOfTicks);			
		}
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		String action = e.getActionCommand();
		if (action.equals(COLLECTION_ACTION)){
			//triggered on collection dropdown action
			JComboBox source = (JComboBox)e.getSource();
			Object selected = source.getSelectedItem();
			if (selected instanceof CyProject){
				CyProject entry = (CyProject)selected;
				model.setSelectedProject(entry);
				this.updateVizButton.setEnabled(true);
			} else {
				model.setSelectedProject(null);
				this.updateVizButton.setEnabled(false);
			}
		} else if (action.equals(MODE_ACTION)){
			JComboBox source = (JComboBox)e.getSource();
			ComparisonMode mode = (Model.ComparisonMode)source.getSelectedItem();
			model.setMode(mode);
		} else if (action.equals(GENERATE_DIFF_ACTION)){
			JCheckBox cb = (JCheckBox)e.getSource();
			model.setGenerateDiffNets(cb.isSelected());
		} else if (action.equals(GENERATE_CONSENSUS_ACTION)){
			JCheckBox cb = (JCheckBox)e.getSource();
			model.setGenerateConsensusNets(cb.isSelected());
		} else if (action.equals(REQUIRE_REF_NET_ACTION)){
			JCheckBox refCb = (JCheckBox)e.getSource();
			model.setRefIncludedInOverlapSupportCutoff(refCb.isSelected());
		}
	}

	@Override
	public void tableChanged(TableModelEvent e) {
		// triggered when the data of the network selection table has changed
		Object source = e.getSource();
		if (source instanceof SelectionTableModel){
			this.refreshCyProject();
			this.updateSupportSlider();
			this.filterTableModel.refresh(model.getSelectedProject());
		} else if (source instanceof EdgeFilterTableModel){
		}
		applyEdgeFilters();			
	}
	
	/**
	 * Read the data from the {@link EdgeFilterTableModel} and apply those filters on 
	 * all Views in the current project.
	 * 
	 */
	private void applyEdgeFilters() {
		Set<String> hiddenInteractions = filterTableModel.getHiddenInteractions();
		Set<CyNetworkView> allViews = new HashSet<CyNetworkView>();
		allViews.addAll(model.getSelectedProject().getAllSourceViews(model.getServices()));
		allViews.addAll(model.getSelectedProject().getAllDifferentialViews(model.getServices()));
		
		System.out.println("-----");
		for (String intt : hiddenInteractions){
			System.out.println(intt);
		}
		System.out.println("-----");
		
		
		for (CyNetworkView cyView : allViews){
			CyNetwork cyNetwork = cyView.getModel();
			CyTable edgeTable = cyNetwork.getDefaultEdgeTable();
			
			for (View<CyNode> nodeView : cyView.getNodeViews()){
				CyNode node = nodeView.getModel();
				Iterable<CyEdge> edgesIt = cyNetwork.getAdjacentEdgeIterable(node, Type.ANY);
				int visibleEdgeCount = 0;
				for (CyEdge edge : edgesIt){
					Long edgeSUID = edge.getSUID();
					CyRow edgeRow = edgeTable.getRow(edgeSUID);
					String interaction = edgeRow.get("interaction", String.class);
					boolean isHidden = hiddenInteractions.contains(interaction);
					View<CyEdge> edgeView = cyView.getEdgeView(edge);
					edgeView.setVisualProperty(BasicVisualLexicon.EDGE_VISIBLE, !isHidden);
					if (!isHidden){
						visibleEdgeCount++;
					}
				}
				nodeView.setVisualProperty(BasicVisualLexicon.NODE_VISIBLE, visibleEdgeCount!=0);
			}
			
			cyView.updateView();
		
		}
	}

	/**
	 * Reads the selections in the GUI and reflects them in the current {@link CyProject} object.
	 */
	public void refreshCyProject(){
		CyProject cyProject = model.getSelectedProject();
		if (cyProject !=null){
			cyProject.setReferenceNetwork(this.networkTableModel.getReferenceNetwork());
			cyProject.setConditionalNetworks(this.networkTableModel.getConditionalNetworks());
			this.runButton.setEnabled(cyProject.canExecute(model));
			
			UpdateVisualStyleTaskFactory tf = new UpdateVisualStyleTaskFactory(model, cyProject);
			TaskIterator it = tf.createTaskIterator();
			DialogTaskManager dtm = model.getServices().getDialogTaskManager();
			dtm.execute(it);					
		}
	}

	@Override
	public void stateChanged(ChangeEvent e) {
		Object source = e.getSource();
		if (source instanceof JSpinner){
			//triggered when cutoff spinner is changed
			JSpinner spinner = (JSpinner)source;
			this.model.setCutoff((Double)spinner.getValue());			
		} else if (source instanceof JSlider){
			JSlider slider = (JSlider)source;
			if (!slider.getValueIsAdjusting()){
				this.model.setOverlapSupportCutoff(slider.getValue());				
			}
		}
		
	}

	@Override
	public void valueChanged(ListSelectionEvent e) {
		//triggered when a row gets selected
		int row = networkTable.getSelectedRow();
		if (row>=0){
			NetworkEntry selectedEntry = networkTableModel.getNetworkEntry(row);
			CyNetwork selectedNetwork = selectedEntry.getNetwork();
			CyNetworkViewManager viewManager = model.getServices().getCyNetworkViewManager();
			Collection<CyNetworkView> cyViews = viewManager.getNetworkViews(selectedNetwork);
			for (CyNetworkView cyView : cyViews){
				model.getServices().getCyApplicationManager().setCurrentNetworkView(cyView);
			}
			networkTable.setRowSelectionInterval(row, row);
		}
		
	}

	
	
}
