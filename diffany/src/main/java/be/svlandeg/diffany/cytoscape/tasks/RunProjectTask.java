package be.svlandeg.diffany.cytoscape.tasks;

import java.util.Collection;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.subnetwork.CyRootNetwork;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskMonitor;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
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
		
		this.runAlgorithm();
		
		taskMonitor.setProgress(1.0);
	}
	
	
	
	private void runAlgorithm() throws InvalidProjectException{
		CyProject cyProject = model.getCurrentProject();
		Project project = cyProject.getProject();
		
		//TODO get cutoff from gui slider
		double cutoff = 0.25;
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(project, cutoff);
		
		addDifferentialNetworks(project.getDifferentialNetworks(), cyProject);
		
	}

	
	private CyNetworkView addCyNetwork(CyNetworkBridge bridge, Network network){
		CyRootNetwork collection = model.getSelectedCollection();
		
		CyNetwork cyNet = bridge.createCyNetwork(network, collection);
		model.getServices().getCyNetworkManager().addNetwork(cyNet);
		
		CyNetworkView cyView = model.getServices().getCyNetworkViewFactory().createNetworkView(cyNet);
		model.getServices().getCyNetworkViewManager().addNetworkView(cyView);
		
		return cyView;
	}
	
	private void addDifferentialNetworks(Collection<DifferentialNetwork> differentialNetworks, CyProject cyProject) {
		CyNetworkBridge bridge = new CyNetworkBridge();
		for (DifferentialNetwork network : differentialNetworks){
			//add the diffnet
			CyNetworkView cyDiffView = this.addCyNetwork(bridge, network);
			
			model.getDiffStyle().apply(cyDiffView);
			cyDiffView.updateView();
			
			//add the overlap
			OverlappingNetwork overlap = network.getOverlappingNetwork();
			CyNetworkView cyOverlapView = this.addCyNetwork(bridge, overlap);
			
			
			model.getSourceStyle().apply(cyOverlapView);
			cyOverlapView.updateView();
			
			cyProject.addResultPair(cyDiffView.getModel(), cyOverlapView.getModel());
		}
		
	}
}
