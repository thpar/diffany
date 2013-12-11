package be.svlandeg.diffany.cytoscape.tasks;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.work.Task;
import org.cytoscape.work.TaskMonitor;

import be.svlandeg.diffany.algorithms.CalculateDiff;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.OverlappingNetwork;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.cytoscape.CyNetworkBridge;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.DefaultNodeMapper;

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
		
		this.runAlgorithm();
		
		taskMonitor.setProgress(1.0);
	}
	
	private void runAlgorithm(){
		CyNetworkBridge bridge = new CyNetworkBridge(model);
		CyNetwork cyRefNetwork = model.getGuiModel().getReferenceNetwork();
		ReferenceNetwork refNet = bridge.getReferenceNetwork(cyRefNetwork);
		
		System.out.println("Refnet");
		System.out.println(refNet.getStringRepresentation());
		System.out.println(refNet.writeEdgesTab());
		
		Set<CyNetwork> cyCondNetworks = model.getGuiModel().getConditionEntries();
		
		System.out.println("Condnets");
		Set<ConditionNetwork> condSet = new HashSet<ConditionNetwork>();
		for (CyNetwork cyNet : cyCondNetworks){
			ConditionNetwork condNet = bridge.getConditionNetwork(cyNet);
			condSet.add(condNet);			
			System.out.println(condNet.getStringRepresentation());
			System.out.println(condNet.writeEdgesTab());
			System.out.println("---");
		}
		
		
		
		Project newProject = new Project("Default Project", refNet, condSet, new DefaultEdgeOntology(), new DefaultNodeMapper());
		
		double cutoff = 0.25;
		new CalculateDiff().calculateAllPairwiseDifferentialNetworks(newProject, cutoff);
		
		addDifferentialNetworks(newProject.getDifferentialNetworks());
	}

	
	private void addDifferentialNetworks(Collection<DifferentialNetwork> differentialNetworks) {
		VisualStyle sourceStyle = model.getGuiModel().getSourceStyle();
		VisualStyle diffStyle = model.getGuiModel().getDiffStyle();
		
		System.out.println("Diffnets");
		CyNetworkBridge bridge = new CyNetworkBridge(model);
		for (DifferentialNetwork network : differentialNetworks){
			System.out.println(network.getStringRepresentation());
			System.out.println(network.writeEdgesTab());
			
			//add the diffnet
			CyNetworkView cyDiffView = this.addCyNetwork(bridge, network);
			sourceStyle.apply(cyDiffView);
			cyDiffView.updateView();
			
			//add the overlap
			OverlappingNetwork overlap = network.getOverlappingNetwork();
			CyNetworkView cyOverlapView = this.addCyNetwork(bridge, overlap);
			diffStyle.apply(cyOverlapView);
			cyOverlapView.updateView();
			
		}
	}
	
	private CyNetworkView addCyNetwork(CyNetworkBridge bridge, Network network){
		CyNetwork cyNet = bridge.createCyNetwork(network);
		model.getServices().getCyNetworkManager().addNetwork(cyNet);
		CyNetworkView cyView = model.getServices().getCyNetworkViewFactory().createNetworkView(cyNet);
		model.getServices().getCyNetworkViewManager().addNetworkView(cyView);
		return cyView;
	}
}
