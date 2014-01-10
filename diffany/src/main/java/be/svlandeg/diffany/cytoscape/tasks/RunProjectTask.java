package be.svlandeg.diffany.cytoscape.tasks;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collection;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import org.cytoscape.application.swing.CySwingApplication;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.model.subnetwork.CySubNetwork;
import org.cytoscape.view.layout.CyLayoutAlgorithm;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskManager;
import org.cytoscape.work.TaskMonitor;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Logger;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.OverlappingNetwork;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.cytoscape.CyNetworkBridge;
import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.InvalidProjectException;
import be.svlandeg.diffany.cytoscape.Model;

/**
 * This Task gathers information from the model and creates and runs a new {@link Project}
 * 
 * @author Thomas Van Parys
 *
 */
public class RunProjectTask implements Task {

	private Model model;
	private CyProject cyProject;
	
	/**
	 * Constructs a new run task.
	 * 
	 * @param model the Diffany {@link Model}
	 * @param cyProject the {@link CyProject} to run.
	 */
	public RunProjectTask(Model model, CyProject cyProject) {
		this.model = model;
		this.cyProject = cyProject;
	}

	@Override
	public void cancel() {
		System.out.println("Task cancelled");
	}

	@Override
	public void run(TaskMonitor taskMonitor) throws Exception {
		taskMonitor.setTitle("Run Diffany Project");
		taskMonitor.setProgress(0.1);
		
		this.runAlgorithm();
		
		taskMonitor.setProgress(1.0);
		
	}
	
	
	/**
	 * Read the log from the {@link Project} and display it as a dialog.
	 * 
	 * @param logger
	 */
	private void displayReport(Logger logger) {
		final JDialog reportDialog = new JDialog(model.getParentWindow(), "Report", false);
		reportDialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		
		StringBuffer reportContent = new StringBuffer();
		for (String msg : logger.getAllLogMessages()){
			reportContent.append(msg);
			reportContent.append(System.getProperty("line.separator"));
		}
		
		JPanel logPanel = new JPanel();
		logPanel.setLayout(new BorderLayout());
		JTextArea text = new JTextArea(reportContent.toString());
		text.setEditable(false);
		text.setLineWrap(true);
		text.setRows(20);
		text.setColumns(50);
		JScrollPane scrollPane = new JScrollPane(text);
		logPanel.add(scrollPane, BorderLayout.CENTER);
		
		JPanel buttonPanel = new JPanel();
		JButton closeButton = new JButton("Close");
		buttonPanel.add(closeButton);
		closeButton.addActionListener(new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent e) {
				reportDialog.dispose();
			}
		});
		logPanel.add(buttonPanel, BorderLayout.SOUTH);
		
		reportDialog.setContentPane(logPanel);
		reportDialog.pack();
		reportDialog.setVisible(true);
		
	}

	/**
	 * Select the correct run mode and run the {@link CyProject}
	 * 
	 * @throws InvalidProjectException thrown when a {@link Project} can not yet be constructed from 
	 * the information in the {@link CyProject}
	 */
	private void runAlgorithm() throws InvalidProjectException{
		Project project = cyProject.getProject();
		
		switch(cyProject.getMode()){
		case REF_PAIRWISE:
			new CalculateDiff().calculateAllPairwiseDifferentialNetworks(project, cyProject.getCutoff());
			break;
		case REF_TO_ALL:	
			new CalculateDiff().calculateOneDifferentialNetwork(project, cyProject.getCutoff());
			break;
		}
		
		addDifferentialNetworks(project.getDifferentialNetworks(), cyProject);
		displayReport(project.getLogger());
	}

	/**
	 * Converts a {@link Network} to a {@link CyNetwork} and registers it
	 * with Cytoscape as a {@link CySubNetwork} of the currently selected Network Collection.
	 * 
	 * 
	 * @param bridge the {@link Network} to {@link CyNetwork} convertor
	 * @param network the {@link Network} to be added
	 * @return the created and added {@link CyNetwork}
	 */
	private CyNetwork addCyNetwork(CyNetworkBridge bridge, Network network){
		CyRootNetwork collection = model.getSelectedCollection();
		
		CyNetwork cyNet = bridge.createCyNetwork(network, collection);
		model.getServices().getCyNetworkManager().addNetwork(cyNet);
		
		CyNetworkView cyView = model.getServices().getCyNetworkViewFactory().createNetworkView(cyNet);
		model.getServices().getCyNetworkViewManager().addNetworkView(cyView);
		
		CyLayoutAlgorithm layout = model.getServices().getCyLayoutAlgorithmManager().getLayout("force-directed");
		TaskIterator it = layout.createTaskIterator(cyView, layout.createLayoutContext(), CyLayoutAlgorithm.ALL_NODE_VIEWS, null);
		TaskManager tm = model.getServices().getTaskManager();
		tm.execute(it);
		
		return cyNet;
	}
	
	/**
	 * Convert and add the resulting networks after running the algorithm.
	 * 
	 * @param differentialNetworks
	 * @param cyProject
	 */
	private void addDifferentialNetworks(Collection<DifferentialNetwork> differentialNetworks, CyProject cyProject) {
		CyNetworkBridge bridge = new CyNetworkBridge();
		for (DifferentialNetwork network : differentialNetworks){
			
			//add the diffnet
			CyNetwork cyDiffNet = this.addCyNetwork(bridge, network);
						
			//add the overlap
			OverlappingNetwork overlap = network.getOverlappingNetwork();
			CyNetwork cyOverlapNet = this.addCyNetwork(bridge, overlap);
			
			cyProject.addResultPair(cyDiffNet, cyOverlapNet);
		}
		
	}
}
