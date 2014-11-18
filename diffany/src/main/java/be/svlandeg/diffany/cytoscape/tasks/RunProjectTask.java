package be.svlandeg.diffany.cytoscape.tasks;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import org.cytoscape.work.Task;
import org.cytoscape.work.TaskMonitor;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;
import be.svlandeg.diffany.core.listeners.ExecutionProgress;
import be.svlandeg.diffany.core.project.LogEntry;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.project.RunConfiguration;
import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.InvalidRunConfigurationException;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.Model.OverlapOperator;

/**
 * This Task gathers information from the model and runs the {@link Project}
 * 
 * @author Thomas Van Parys
 *
 */
public class RunProjectTask extends ExecutionProgress implements Task  {

	private Model model;
	private CyProject cyProject;
	private TaskMonitor taskMonitor;
	
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
		this.taskMonitor = taskMonitor;
		this.taskMonitor.setTitle("Diffany Algorithm");
		this.taskMonitor.setProgress(0.01);
		
		this.runAlgorithm();
		
		this.taskMonitor.setProgress(1.0);
		
	}
	
	/**
	 * Generates the {@link RunConfiguration}, selects the correct run mode and runs the {@link Project}
	 * 
	 * @throws InvalidRunConfigurationException thrown when a {@link RunConfiguration} can not yet be constructed from 
	 * the information in the {@link CyProject}
	 */
	private void runAlgorithm() throws InvalidRunConfigurationException{
		int runId = cyProject.generateRunConfiguration(model, this);
		
		int nextID = 666;
		int generateDiff = -1;
		if (model.isGenerateDiffNets())
		{
			generateDiff = nextID;
			nextID++;
		}
		int generateConsensus = -1;
		if (model.isGenerateConsensusNets())
		{
			generateConsensus = nextID;
			nextID++;
		}
		
		
		switch(model.getMode()){
		case REF_PAIRWISE:
			new CalculateDiff().calculateAllPairwiseDifferentialNetworks(cyProject.getProject(), runId, 
					model.getCutoff(), model.isGenerateDiffNets(), model.isGenerateConsensusNets(), 
					nextID, model.getOverlapOperator() == OverlapOperator.MIN, this);
			break;
		case REF_TO_ALL:	
			new CalculateDiff().calculateOneDifferentialNetwork(cyProject.getProject(), runId, 
					model.getCutoff(), generateDiff, generateConsensus, model.getOverlapOperator() == OverlapOperator.MIN, this);
			break;
		}
		
		cyProject.update(model.getServices());
		displayReport(cyProject.getProject().getLogger(runId));
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
		for (LogEntry msg : logger.getAllLogMessages()){
			reportContent.append(msg.toString());
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

	@Override
	public void setProgress(String message, int progress, int total) {
		this.taskMonitor.setStatusMessage(message);
		this.taskMonitor.setProgress((double)progress/(double)total);
	}


	
	
}
