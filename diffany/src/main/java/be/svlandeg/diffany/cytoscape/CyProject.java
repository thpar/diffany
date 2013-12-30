package be.svlandeg.diffany.cytoscape;

import java.util.HashSet;
import java.util.Set;

import org.cytoscape.model.CyNetwork;

import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.Project;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualDiffStyle;
import be.svlandeg.diffany.cytoscape.vizmapper.VisualSourceStyle;
import be.svlandeg.diffany.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.semantics.DefaultNodeMapper;



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
	private Set<CyNetwork> conditionalNetworks = new HashSet<CyNetwork>();
	private Set<CyNetwork> resultNetworks = new HashSet<CyNetwork>();;

	
	private VisualSourceStyle visualSourceStyle;
	private VisualDiffStyle visualDiffStyle;



	public CyProject(){
	}
	
	
	
	public CyNetwork getReferenceNetwork() {
		return referenceNetwork;
	}



	public void setReferenceNetwork(CyNetwork referenceNetwork) {
		this.referenceNetwork = referenceNetwork;
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
		
		CyNetworkBridge bridge = new CyNetworkBridge();
		ReferenceNetwork refNet = bridge.getReferenceNetwork(this.getReferenceNetwork());
		
		Set<ConditionNetwork> condNets = new HashSet<ConditionNetwork>();
		for (CyNetwork cyCondNet : this.getConditionalNetworks()){
			ConditionNetwork condNet = bridge.getConditionNetwork(cyCondNet);
			condNets.add(condNet);			
		}
		
		this.project = new Project("New Project", refNet, condNets, new DefaultEdgeOntology(), new DefaultNodeMapper());
		
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
