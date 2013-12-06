package be.svlandeg.diffany.cytoscape.tasks;

import org.cytoscape.work.Task;
import org.cytoscape.work.TaskMonitor;

import be.svlandeg.diffany.cytoscape.Model;

/**
 * Temporary task to experiment.
 * 
 * @author thpar
 *
 */
public class RunProjectTask implements Task {

	private Model model;
	
	
	public RunProjectTask(Model model) {
		this.model = model;
	}

	@Override
	public void cancel() {
		System.out.println("Task cancelled");
	}

	@Override
	public void run(TaskMonitor taskMonitor) throws Exception {
		taskMonitor.setTitle("Run Diffany Project");
		taskMonitor.setProgress(0.1);
		
		model.runAlgorithm();
		
		taskMonitor.setProgress(1.0);
	}

}
