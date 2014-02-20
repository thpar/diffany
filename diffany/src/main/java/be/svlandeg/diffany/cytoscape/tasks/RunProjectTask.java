package be.svlandeg.diffany.cytoscape.tasks;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collection;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

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
import be.svlandeg.diffany.cytoscape.InvalidRunConfigurationException;
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
	 * @throws InvalidRunConfigurationException thrown when a {@link Project} can not yet be constructed from 
	 * the information in the {@link CyProject}
	 */
	private void runAlgorithm() throws InvalidRunConfigurationException{
		int runId = cyProject.generateRunConfiguration();
		
		switch(model.getMode()){
		case REF_PAIRWISE:
			new CalculateDiff().calculateAllPairwiseDifferentialNetworks(cyProject.getProject(), runId, model.getCutoff());
			break;
		case REF_TO_ALL:	
			new CalculateDiff().calculateOneDifferentialNetwork(cyProject.getProject(), runId, model.getCutoff());
			break;
		}
		
		cyProject.update(model.getServices());
		displayReport(cyProject.getProject().getLogger(runId));
	}

	
	
}
