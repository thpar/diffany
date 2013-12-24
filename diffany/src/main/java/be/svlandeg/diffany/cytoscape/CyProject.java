package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.Set;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.view.model.CyNetworkView;

import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualDiffStyle;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualSourceStyle;



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
	private Set<CyNetwork> resultNetworks;

	
	private VisualSourceStyle visualSourceStyle;
	private VisualDiffStyle visualDiffStyle;
	

	private Model model;

	public CyProject(Model model){
		this.model = model;
	}
	
	
	
	public CyNetwork getReferenceNetwork() {
		return referenceNetwork;
	}



	public void setReferenceNetwork(CyNetwork referenceNetwork) {
		this.referenceNetwork = referenceNetwork;
		Collection<CyNetworkView> networkViews = model.getServices().getCyNetworkViewManager().getNetworkViews(referenceNetwork);
		for (CyNetworkView view : networkViews){
			visualSourceStyle.getVisualStyle().apply(view);
		}
	}

	public void removeConditionalNetwork(CyNetwork cyNetwork){
		this.conditionalNetworks.remove(cyNetwork);
	}


	public Set<CyNetwork> getConditionalNetworks() {
		return conditionalNetworks;
	}



	public void setConditionalNetworks(Set<CyNetwork> conditionalNetworks) {
		this.conditionalNetworks = conditionalNetworks;
	}

	public void addConditionalNetwork(CyNetwork conditionalNetwork){
		this.conditionalNetworks.add(conditionalNetwork);
		Collection<CyNetworkView> networkViews = model.getServices().getCyNetworkViewManager().getNetworkViews(conditionalNetwork);
		for (CyNetworkView view : networkViews){
			visualSourceStyle.getVisualStyle().apply(view);
		}
	}


	public Set<CyNetwork> getDifferentialNetworks() {
		return resultNetworks;
	}



	public void setDifferentialNetworks(Set<CyNetwork> resultNetworks) {
		this.resultNetworks = resultNetworks;
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


	
	/**
	 * Gets the visual style for source (input and overlap) networks, wrapped in a {@link VisualSourceStyle}
	 * and adjusted to the edge types used in the networks of the current {@link CyProject}
	 * @return the visual style for source networks
	 */
	public VisualSourceStyle getVisualSourceStyle() {
		return visualSourceStyle;
	}


	/**
	 * Gets the visual style for differential networks, wrapped in a {@link VisualDiffStyle}
	 * and adjusted to the edge types used in the networks of the current {@link CyProject}
	 * @return the visual style for differential networks
	 */
	public VisualDiffStyle getVisualDiffStyle() {
		return visualDiffStyle;
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
	}



	public void addResultNetwork(CyNetwork cyNet) {
		this.resultNetworks.add(cyNet);
		
	}
		

}
