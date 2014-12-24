package be.svlandeg.diffany.core.algorithms;

/*
 * #%L
 * Diffany
 * %%
 * Copyright (C) 2014 PSB/UGent - Sofie Van Landeghem and Thomas Van Parys
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */


import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import be.svlandeg.diffany.core.networks.Node;

/**
 * This class defines the nodes to be used in output (differential/consensus) networks, by transferring the node attributes and name/ID from the corresponding input nodes.
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
	 * @throws IllegalArgumentException when both the refnode is null and the set of conditional nodes is also null or empty
	 */
	public Node getConsensusNode(Node refNode, Set<Node> conNodes, int supportingCutoff) throws IllegalArgumentException
	{
		if (refNode == null && (conNodes == null || conNodes.isEmpty()))
		{
			String errormsg = "At least one non-null node should be specified! ";
			throw new IllegalArgumentException(errormsg);
		}
		
		Set<Node> allNodes = new HashSet<Node>();
		if (conNodes != null)
		{
			allNodes.addAll(conNodes);
		}

		Node exampleNode = refNode;
		Iterator<Node> it = null;
		if (refNode == null)
		{
			it = conNodes.iterator();
			exampleNode = it.next();
		}
		else
		{
			allNodes.add(refNode);
		}
		Node consensusNode = new Node(exampleNode.getID(), exampleNode.getDisplayName());

		// collect all known node attributes from the input data
		Set<String> attributes = new HashSet<String>();
		for (Node n : allNodes)
		{
			attributes.addAll(n.getAllAttributeNames());
		}

		// calculate DE status from the condition-specific nodes alone
		String de_attribute = Node.de_attribute;
		if (attributes.contains(de_attribute))
		{
			attributes.remove(de_attribute);
			int countUp = 0;
			int countDown = 0;

			for (Node n : conNodes)
			{
				String deStatus = (String) n.getAttribute(de_attribute);
				if (deStatus != null && deStatus.equals(Node.upregulated))
				{
					countUp++;
				}
				if (deStatus != null && deStatus.equals(Node.downregulated))
				{
					countDown++;
				}
			}
			if (countUp > countDown && countUp >= supportingCutoff)
			{
				consensusNode.setAttribute(de_attribute, Node.upregulated);
			}
			else if (countDown > countUp && countDown >= supportingCutoff)
			{
				consensusNode.setAttribute(de_attribute, Node.downregulated);
			}
			else
			{
				// We also put it to 'no' when the count of ups and downs is equal, because then there is no consensus
				consensusNode.setAttribute(de_attribute, Node.not_de);
			}
		}

		// for each (non-DE) attribute, keep searching until a non-null value is found (if ever)
		for (String att_name : attributes)
		{
			Object att_value = exampleNode.getAttribute(att_name);

			while (att_value == null && it != null && it.hasNext())
			{
				att_value = it.next().getAttribute(att_name);
			}
			if (att_value != null)
			{
				consensusNode.setAttribute(att_name, att_value);
			}
		}

		return consensusNode;
	}

}
