package be.svlandeg.diffany.cytoscape;

import java.util.Collection;

import org.cytoscape.model.CyNetwork;
import org.cytoscape.model.CyNetworkFactory;
import org.cytoscape.model.CyNetworkManager;
import org.cytoscape.model.CyNode;
import org.cytoscape.model.CyRow;
import org.cytoscape.model.CyTable;

import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Network;
import be.svlandeg.diffany.concepts.Node;

/**
 * Conversion class that takes the appropriate Cytoscape Factories and uses them
 * to create given concepts.
 * 
 * @author thpar
 *
 */
public class CyNetworkBridge {
	
	private CyNetworkFactory cyNetworkFactory;
	
	private final String NORMALIZED_NAME = "normalized_name";

	public CyNetworkBridge(CyNetworkFactory cyNetworkFactory){
		this.cyNetworkFactory = cyNetworkFactory;
	}

	/**
	 * Converts a {@link Network} to a {@link CyNetwork}
	 * 
	 * @param network input {@link Network} 
	 * @return {@link CyNetwork} that can be added to Cytoscape using a {@link CyNetworkManager}
	 */
	public CyNetwork createCyNetwork(Network network){
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
	 * 
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
}
