package be.svlandeg.diffany.core.algorithms;

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
	 * Method that defines the differential node from the corresponding nodes in the reference and condition-specific networks.
	 * The IDs of all the nodes involved should be the same, as well as the name. Attributes are merged.
	 * 
	 * The main functionality of this method is then to calculate the DE status of the differential node.
	 * 
	 * An important parameter is supportingCutoff, which determines the amount of support needed for a node to be considered as differentially expressed (DE). 
	 * If it equals the number of input networks, all nodes need to agree on the DE status.
	 * However, if it is smaller, e.g. 3 out of 4, there can be one 'outlier' node, allowing some noise in the input and creating more robust interpretations.
	 * The supportingCutoff should ideally be somewhere between 50% and 100%, but this choice is determined by the specific use-case / application. Instead of being a percentage, this method requires the support to be expressed
	 * as a minimal number of supporting nodes (networks).
	 * 
	 * @param conNodes the nodes in the condition-specific networks
	 * @param supportingCutoff the minimal number of condition networks (inclusive) that need to reach the same DE consensus for the differential node to be DE as well.
	 * 
	 * @return the node in the differential network, with its proper attributes.
	 * @throws IllegalArgumentException when the type of the reference or condition-specific edge does not exist in this ontology
	 */
	public Node getConsensusNode(Set<Node> conNodes, int supportingCutoff) throws IllegalArgumentException
	{
		// TODO
		return null;
		
		// TODO: ref required?
	}
	

}
