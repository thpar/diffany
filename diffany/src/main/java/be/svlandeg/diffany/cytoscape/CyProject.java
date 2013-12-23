package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.Set;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.vizmap.VisualStyle;

import be.svlandeg.diffany.concepts.DifferentialNetwork;
import be.svlandeg.diffany.concepts.OverlappingNetwork;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualDiffStyle;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualSourceStyle;
import be.svlandeg.diffany.internal.Services;



/**
 * The CyProject keeps the GUI settings for the current project and is able to generate a 
 * {@link Project} object to execute the algorithm.
 * 
 * @author thpar
 *
 */
public class CyProject{

	Project project;
	
	private CyNetwork referenceNetwork;
	private Set<CyNetwork> conditionalNetworks;
	private Set<CyNetwork> differentialNetworks;

	
	private VisualSourceStyle visualSourceStyle;
	private VisualDiffStyle visualDiffStyle;
	
	
	private Services services;

	private Model model;

	public CyProject(Model model){
		this.services = model.getServices();
		this.model = model;
	}
	
	
	
	public CyNetwork getReferenceNetwork() {
		return referenceNetwork;
	}



	public void setReferenceNetwork(CyNetwork referenceNetwork) {
		this.referenceNetwork = referenceNetwork;
	}



	public Set<CyNetwork> getConditionalNetworks() {
		return conditionalNetworks;
	}



	public void setConditionalNetworks(Set<CyNetwork> conditionalNetworks) {
		this.conditionalNetworks = conditionalNetworks;
	}

	public void addConditionalNetwork(CyNetwork conditionalNetwork){
		this.conditionalNetworks.add(conditionalNetwork);
		Collection<CyNetworkView> networkViews = services.getCyNetworkViewManager().getNetworkViews(conditionalNetwork);
		for (CyNetworkView view : networkViews){
			visualSourceStyle.getVisualStyle().apply(view);
		}
	}


	public Set<CyNetwork> getDifferentialNetworks() {
		return differentialNetworks;
	}



	public void setDifferentialNetworks(Set<CyNetwork> differentialNetworks) {
		this.differentialNetworks = differentialNetworks;
	}


	/**
	 * Generates a {@link Project} from current information.
	 * 
	 * @return a newly constructed {@link Project}
	 * @throws InvalidProjectException is thrown when not all necessary parameters are there to construct the {@link Project}
	 */
	public Project getProject() throws InvalidProjectException{
		if (!canExecute()){
			throw new InvalidProjectException();
		}
		
		return project;
	}


	

	public VisualSourceStyle getVisualSourceStyle() {
		return visualSourceStyle;
	}



	public VisualDiffStyle getVisualDiffStyle() {
		return visualDiffStyle;
	}

	private void addDifferentialNetworks(Collection<DifferentialNetwork> differentialNetworks) {
		VisualStyle sourceStyle = model.getGuiModel().getSourceStyle().getVisualStyle();
		VisualStyle diffStyle = model.getGuiModel().getDiffStyle().getVisualStyle();
		
		System.out.println("Diffnets");
		CyNetworkBridge bridge = new CyNetworkBridge(model);
		for (DifferentialNetwork network : differentialNetworks){
			System.out.println(network.getStringRepresentation());
			System.out.println(network.writeEdgesTab());
			
			//add the diffnet
			CyNetworkView cyDiffView = this.addCyNetwork(bridge, network);
			diffStyle.apply(cyDiffView);
			cyDiffView.updateView();
			
			//add the overlap
			OverlappingNetwork overlap = network.getOverlappingNetwork();
			CyNetworkView cyOverlapView = this.addCyNetwork(bridge, overlap);
			sourceStyle.apply(cyOverlapView);
			cyOverlapView.updateView();
			
		}
		
	}

	/**
	 * Checks if all necessary parameters are set to execute the algorithm.
	 * 
	 * @return true when minimal input parameters are set
	 */
	public boolean canExecute(){
		//TODO actual check 
		return true;
	}


	/**
	 * Iterates all resulting networks (differential and overlap) and creates {@link CyNetwork}s and {@link CyView}s
	 * where needed.
	 */
	public void updateResultNetworks() {
		// TODO Auto-generated method stub
		
	}
		

}
