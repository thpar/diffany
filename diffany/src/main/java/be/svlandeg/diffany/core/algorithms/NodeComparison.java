package be.svlandeg.diffany.core.algorithms;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import be.svlandeg.diffany.core.networks.Node;

/**
 * This class defines the nodes to be used in output (differential/consensus) networks
 * 
 * @author Sofie Van Landeghem
 */
public class NodeComparison
{

	/**
	 * Method that defines the consensus node from the corresponding nodes in the reference and condition-specific networks.
	 * The IDs of all the nodes involved should be the same, as well as the name. Attributes are merged.
	 * 
	 * The main functionality of this method is then to calculate the DE status of the differential node, by combining the information from the condition-specific networks. 
	 * For this, the supportingCutoff is an important parameter and determines the amount of support needed for a node to be considered as differentially expressed (DE). 
	 * If it equals the number of input networks, all nodes need to agree on the DE status (except for the reference node, which is non-DE by default).
	 * 
	 * However, if the cutoff is smaller, e.g. 3 out of 4, there can be one 'outlier' node, allowing some noise in the input and creating more robust interpretations.
	 * The supportingCutoff should ideally be somewhere between 50% and 100%, but this choice is determined by the specific use-case / application. Instead of being a percentage, this method requires the support to be expressed
	 * as a minimal number of supporting nodes (networks).
	 * 
	 * @param refNode the node in the reference network - may be null in case no reference is specified/known
	 * @param conNodes the nodes in the condition-specific networks
	 * @param supportingCutoff the minimal number of condition networks (inclusive) that need to reach the same DE consensus for the differential node to be DE as well
	 * 
	 * @return the node in the differential network, with its proper attributes.
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public Node getConsensusNode(Node refNode, Set<Node> conNodes, int supportingCutoff) throws IllegalArgumentException
	{
		Set<Node> allNodes = new HashSet<Node>();
		allNodes.addAll(conNodes);
		
		Node exampleNode = refNode;
		Iterator<Node> it = conNodes.iterator();
		if (refNode == null)
		{
			exampleNode = it.next();
		}
		else
		{
			allNodes.add(refNode);
		}
		Node consensusNode = new Node(exampleNode.getID(), exampleNode.getDisplayName(false));
		
		
		// for all attributes except the DE property, it does not matter which is reference and which condition-specific
		Set<String> attributes = new HashSet<String>();
		for (Node n : allNodes)
		{
			attributes.addAll(n.getAllAttributeNames());
		}
		String de_attribute = Node.de_attribute;
		attributes.remove(de_attribute);
		
		// for each (non-DE) attribute, keep searching until a non-null value is found (if ever)
		for (String att_name : attributes)
		{
			String att_value = exampleNode.getAttribute(att_name);
			
			while (att_value == null && it.hasNext())
			{
				att_value = it.next().getAttribute(att_name);
			}
			if (att_value != null)
			{
				consensusNode.setAttribute(att_name, att_value);
			}
		}
		
		// TODO: calculate DE status
		
		
		return consensusNode;
	}
	

}
