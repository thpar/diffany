package be.svlandeg.diffany.cytoscape;

import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.cytoscape.model.CyEdge;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyRow;
import org.cytoscape.model.CyTable;

import be.svlandeg.diffany.concepts.Condition;
import be.svlandeg.diffany.concepts.ConditionNetwork;
import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;
import be.svlandeg.diffany.concepts.ReferenceNetwork;
import be.svlandeg.diffany.internal.Services;

/**
 * Conversion class that takes the appropriate Cytoscape Factories and uses them
 * to create given concepts. The other way around, it enables existing {@link CyNetwork}s to
 * be converted to the project's {@link Network} concept. 
 * 
 * @author thpar
 *
 */
public class CyNetworkBridge {
	

	/**
	 * Column name of an extra field to be used to map nodes from different networks onto each other.
	 */
	private final String NORMALIZED_NAME = "normalized_name";

	private Services services;

	private Model model;

	public CyNetworkBridge(Model model){
		this.services = model.getServices();
		this.model = model;
	}

	/**
	 * Converts a {@link Network} to a {@link CyNetwork}
	 * 
	 * @param network input {@link Network} 
	 * @return {@link CyNetwork} that can be added to Cytoscape using a {@link CyNetworkManager}
	 */
	public CyNetwork createCyNetwork(Network network){
		CyNetworkFactory cyNetworkFactory = services.getCyNetworkFactory();
		CyNetwork cyNetwork = cyNetworkFactory.createNetwork();
		cyNetwork.getRow(cyNetwork).set(CyNetwork.NAME,network.getName());
		cyNetwork.getDefaultNodeTable().createColumn(NORMALIZED_NAME, String.class, false);
		
		
		for (Node node: network.getNodes()){
			CyNode cyNode = cyNetwork.addNode();
			cyNetwork.getRow(cyNode).set(CyNetwork.NAME, node.getName());
			cyNetwork.getRow(cyNode).set(NORMALIZED_NAME, node.getName(true));
		}
		
		
		for (Edge edge : network.getEdges()){
			Node from = edge.getSource();
			Node to = edge.getTarget();
			
			CyTable table = cyNetwork.getDefaultNodeTable();
			final String pk = table.getPrimaryKey().getName();
			CyNode fromNode = nodeToCyNode(from, cyNetwork, table, pk);
			CyNode toNode = nodeToCyNode(to, cyNetwork, table, pk);
			
			cyNetwork.addEdge(fromNode, toNode, !edge.isSymmetrical());
		}
		return cyNetwork;
	}
	
	/**
	 * Takes a {@link Node} and converts it to a {@link CyNode}, using the normalized_name column.
	 * 
	 * @param node
	 * @param cyNetwork
	 * @param table
	 * @param pk
	 * @return
	 */
	private CyNode nodeToCyNode(Node node, CyNetwork cyNetwork, CyTable table, String pk){
		CyNode cyNode = null;
		Collection<CyRow> matches = table.getMatchingRows(NORMALIZED_NAME, node.getName(true));
		for (final CyRow row : matches){
			Long suid = row.get(pk, Long.class);
			cyNode = cyNetwork.getNode(suid);
		}
		return cyNode;
	}
	
	enum NetworkType{
		REFERENCE, CONDITION;
	}
	
	/**
	 * Gets the {@link CyNetwork} from Cytoscape, based on the shared name,
	 * and creates a {@link Network} out of it.
	 * 
	 * @param networkID the "shared name" of the network in the network table.
	 * 
	 * 
	 * @return the equivalent {@link Network} object
	 */
	public Network getReferenceNetwork(String networkID){
		return getNetwork(networkID, NetworkType.REFERENCE);
	}
	
	/**
	 * Gets the {@link CyNetwork} from Cytoscape, based on the shared name,
	 * and creates a {@link Network} out of it.
	 * 
	 * @param networkID the "shared name" of the network in the network table.
	 * 
	 * 
	 * @return the equivalent {@link Network} object
	 */
	public Network getConditionNetwork(String networkID){
		return getNetwork(networkID, NetworkType.CONDITION);
	}
	
	/**
	 * Gets the {@link CyNetwork} from Cytoscape, based on the name,
	 * and creates a {@link Network} out of it.
	 * 
	 * @param networkID the "name" of the network in the network table.
	 * 
	 * 
	 * @return the equivalent {@link Network} object
	 */
	private Network getNetwork(String networkID, NetworkType type){
		
		Network network = null;
		switch(type){
		case REFERENCE: 
			network = new ReferenceNetwork(networkID);
			break;
		case CONDITION:
			//TODO get conditions from gui
			Set<Condition> conditions = new HashSet<Condition>();
			conditions.add(new Condition("temp_condition"));
			network = new ConditionNetwork(networkID, conditions);
			break;
		}
		
		//get the CyNetwork
		CyNetwork cyNetwork = model.getNetworkByName(networkID);
		List<CyNode> cyNodes = cyNetwork.getNodeList();
		Set<Node> nodeSet = new HashSet<Node>();	
		for (CyNode cyNode : cyNodes){
			
		}
		
		Set<Edge> edgeSet = new HashSet<Edge>();
		List<CyEdge> cyEdges = cyNetwork.getEdgeList();				
		
		
		return network;
	}
	
}
