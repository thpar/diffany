package be.svlandeg.diffany.cytoscape.tasks;

import java.util.Collection;
import java.util.Set;

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
import be.svlandeg.diffany.concepts.ConditionNetwork;
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
	private CyProject cyProject;
	
	
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
	
	
	
	private void runAlgorithm() throws InvalidProjectException{
		Project project = cyProject.getProject();
		
		System.out.println("Running algorithm: "+cyProject.getMode());
		
		switch(cyProject.getMode()){
		case REF_PAIRWISE:
			new CalculateDiff().calculateAllPairwiseDifferentialNetworks(project);
			break;
		case REF_TO_ALL:
			new CalculateDiff().calculateDiffNetwork(project.getReferenceNetwork(), 
					(Set<ConditionNetwork>)project.getConditionNetworks(), 
					project.getEdgeOntology(), project.getNodeMapper());
			break;
		}
		
		addDifferentialNetworks(project.getDifferentialNetworks(), cyProject);
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
