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


import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.ConsensusNetwork;
import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.EdgeDefinition;
import be.svlandeg.diffany.core.networks.InputNetwork;
import be.svlandeg.diffany.core.networks.Network;
import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.progress.ScheduledTask;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.semantics.EdgeOntology;

/**
 * This class provides generic methods useful for network cleaning before or
 * after applying differential algorithms.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkCleaning
{

	private Logger logger;

	/**
	 * Create a new cleaning object, which can log important messages.
	 * 
	 * @param logger the logger object
	 */
	public NetworkCleaning(Logger logger)
	{
		if (logger == null)
		{
			String errormsg = "The logger should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.logger = logger;
	}

	/**
	 * Clean an input network:
	 * Per pair of nodes, group all input edges into subclasses per root category of the EdgeOntology, unify the directionality
	 * (either all symmetric or all directed, as dicated by the edge ontology), and resolve conflicts within a root category.
	 * Resolving conflicts by weight is only an option for input/consensus networks as it uses the source categories of the edge ontology.
	 * 
	 * Be aware: this function changes the input network object!
	 * 
	 * @param net the network that needs cleaning
	 * @param eo the edge ontology
	 * @param task the task object that keeps track of the progress of this calculation (can be null)
	 * @param resolveConflictsByWeight if true, conflicts are resolved by selecting the edge with the highest weight - only possible for Input Networks!
	 * @param unifyEdgeTypes whether or not to unify the edge types - should only be true for input and consensus networks!
	 */
	protected void fullCleaning(Network net, EdgeOntology eo, ScheduledTask task, boolean resolveConflictsByWeight, boolean unifyEdgeTypes)
	{
		boolean isSourceNetwork = net instanceof InputNetwork || net instanceof ConsensusNetwork;
		if (resolveConflictsByWeight && ! isSourceNetwork)
		{
			String errormsg = "Conflicts can not be resolved by weight for this type of network: " + net.getClass();
			throw new IllegalArgumentException(errormsg);
		}
		if (unifyEdgeTypes && ! isSourceNetwork)
		{
			String errormsg = "Edge types can not be unified for this type of network: " + net.getClass();
			throw new IllegalArgumentException(errormsg);
		}
		String progressMessage = "Cleaning network " + net.getName();
		if (task != null)
		{
			task.setMessage(progressMessage);
			task.ticksDone(1);
		}
		
		logger.log(" Full cleaning of " + net.getName());

		// make edges directed when defined as such by the edge ontology
		Set<Node> nodes = net.getNodes();
		Set<Edge> edges = net.getEdges();
		
		if (unifyEdgeTypes)
		{
			edges = new Unification(logger).unifyEdgeDirection(net.getEdges(), eo);
		}
		
		net.setNodesAndEdges(nodes, edges);

		// remove obvious conflicts / redundant edges
		removeRedundantEdges(net, eo, task);
		
		// in case there are still conflicts, keep only the highest weight
		if (resolveConflictsByWeight)
		{
			resolveToOne(net, eo);
		}
		
		if (task != null)
		{
			task.done();
		}
	}

	/**
	 * Clean an output consensus network: Remove redundant/duplicate edges in the network. 
	 * Be aware: this function changes the network object!
	 * 
	 * @param net the consensus network that needs cleaning
	 * @param eo the edge ontology
	 * @param task the task object that keeps track of the progress of this calculation (can be null)
	 */
	public void fullConsensusOutputCleaning(ConsensusNetwork net, EdgeOntology eo, ScheduledTask task)
	{
		fullCleaning(net, eo, task, false, true);
	}

	/**
	 * Clean an output network: Remove redundant/duplicate edges in the network.
	 * 
	 * Be aware: this function changes the network object!
	 * 
	 * @param net the network that needs cleaning
	 * @param eo the edge ontology
	 * @param task the task object that keeps track of the progress of this calculation (can be null)
	 */
	public void fullDifferentialOutputCleaning(Network net, EdgeOntology eo, ScheduledTask task)
	{
		fullCleaning(net, eo, task, false, false);
	}

	/**
	 * Remove redundant edges in the network, such as those that are symmetrical and represented twice (source-target and target-source),
	 * or a generic edge (e.g. regulation) when a specific edge (e.g. inhibition) is also present (of the same weight or more).
	 * 
	 * Otherwise, the weight and negation status should be equal.
	 * 
	 * @param net the network that needs cleaning
	 * @param eo the edge ontology
	 * @param task the task object that keeps track of the progress of this calculation (can be null)
	 */
	protected void removeRedundantEdges(Network net, EdgeOntology eo, ScheduledTask task)
	{
		String progressMessage = "Cleaning network " + net.getName();
		if (task != null)
		{
			task.setMessage(progressMessage);
			task.ticksDone(1);
		}
		logger.log(" Removing redundant edges from network " + net.getName());

		Set<Edge> removed_edges = new HashSet<Edge>();

		Set<Edge> oldEdges = new HashSet<Edge>(net.getEdges());
		int iterations = oldEdges.size();
		int ticksPerReport = 10;
		int iterationPerReport = 0;
		
		if (task != null)
		{
			iterationPerReport = (ticksPerReport * iterations) / task.ticksToGo();
		}
		
		int iterationsDone = 0;
		// remove duplicate symmetrical edges between source-target and target-source
		for (Edge et : oldEdges)
		{
			iterationsDone++;
			if (task != null && iterationPerReport > 0 && iterationsDone % iterationPerReport == 0)
			{
				task.ticksDone(ticksPerReport);
			}
			
			Node n1 = et.getSource();
			Node n2 = et.getTarget();

			Set<Edge> all_edges = net.getAllEdges(n1, n2);
			for (Edge eb : all_edges)
			{
				// comparing two edges 'et' and 'eb', not removed in a previous iteration
				if (!et.equals(eb) && !removed_edges.contains(et) && !removed_edges.contains(eb))
				{
					boolean etParent = false;
					boolean ebParent = false;
					
					String typeET = et.getType();
					String typeEB = eb.getType();
					
					// this check allows the method to be used also for output networks
					if (eo.isDefinedSourceType(typeEB) && eo.isDefinedSourceType(typeET))
					{
						String catEB = eo.getSourceCategory(typeEB);
						String catET = eo.getSourceCategory(typeET);
						
						etParent = eo.isSourceCatChildOf(catEB, catET) > 0;  // if positive, typeET is a parent of typeEB
						ebParent = eo.isSourceCatChildOf(catET, catEB) > 0;  // if positive, typeEB is a parent of typeET 
					}
					
					boolean equalType = typeET.equals(typeEB);
					boolean equalDirection = true;

					// they need to be both symmetrical or both directed
					
					/* Currently this thus does not filter for A -> B which is redundant because of A <-> B
					   but this should never happen because of the edge ontology definitions */
					
					if (et.isSymmetrical() != eb.isSymmetrical())
					{
						equalDirection = false;
					}
					// if they are directed, check that the source & target agree (checking only source should be enough, but just to be sure...)
					if (! et.isSymmetrical())
					{
						if (! et.getSource().getID().equals(eb.getSource().getID()))
						{
							equalDirection = false;
						}
						if (! et.getTarget().getID().equals(eb.getTarget().getID()))
						{
							equalDirection = false;
						}
					}
					
					// both have the same symmetry status and negation status
					if (equalDirection && (et.isNegated() == eb.isNegated()))
					{
						// allow for small rounding errors
						boolean equalWeight = Math.abs(et.getWeight() - Math.abs(eb.getWeight())) < 0.000001;
						
						boolean ebLower = et.getWeight() > eb.getWeight();
						boolean etLower = eb.getWeight() > et.getWeight();
						
						Edge toRemove = null;		
						
						// remove the affirmative parent (or equal) with equal or lower weight
						if ((equalType || ebParent) && ! eb.isNegated() && (equalWeight || ebLower))
						{
							toRemove = eb;
						}
						if ((equalType || etParent) && ! et.isNegated() && (equalWeight || etLower))
						{
							toRemove = et;
						}
						
						// remove the negated child (or equal) if it's of equal weight or lower
						if ((equalType || ebParent) && et.isNegated() && (equalWeight || etLower))
						{
							toRemove = et;
						}
						if ((equalType || etParent) && eb.isNegated() && (equalWeight || ebLower))
						{
							toRemove = eb;
						}
						
						if (toRemove != null)
						{
							net.removeEdge(toRemove);
							removed_edges.add(toRemove);
						}
					}
				}
			}
		}
		if (task != null)
		{
			task.done();
		}
	}

	/**
	 * Create a new list of edges by making either all input interactions directed, or all symmetrical. 
	 * As soon as one input edge is directed, the whole set will become directed.
	 * This method does not impose any other conditions such as whether the edges are between the same nodes or not, of the same type or not.
	 * A sensible grouping of the edges should thus have been done by the calling method.
	 * 
	 * @param oldList the original edges which might have a mixture of directed and symmetrical edges - will not be altered!
	 * @return a new list of edges with their edge directionality unified
	 */
	public List<EdgeDefinition> unifyDirection(List<EdgeDefinition> oldList)
	{
		boolean symmetrical = true;
		for (EdgeDefinition oldEdge : oldList)
		{
			symmetrical = symmetrical && oldEdge.isSymmetrical();
		}
		List<EdgeDefinition> newEdges = new ArrayList<EdgeDefinition>();
		for (EdgeDefinition oldEdge : oldList)
		{
			EdgeDefinition newEdge = new EdgeDefinition(oldEdge);
			newEdge.makeSymmetrical(symmetrical);
			newEdges.add(newEdge);
		}
		return newEdges;
	}


	/**
	 * Clean an input condition-specific network:
	 * Per pair of nodes, group all input edges into subclasses per root category of the EdgeOntology, unify the directionality
	 * (either all symmetric or all directed, as dicated by the edge ontology), and resolve conflicts within a root category.
	 * 
	 * @param net the network that needs cleaning
	 * @param eo the edge ontology
	 * @param task the task object that keeps track of the progress of this calculation (can be null)
	 * 
	 * @return a cleaned condition-specific network representing the same semantic information
	 */
	public ConditionNetwork fullInputConditionCleaning(ConditionNetwork net, EdgeOntology eo, ScheduledTask task)
	{
		ConditionNetwork resultNet = new ConditionNetwork(net.getName(), net.getID(), net.getAllNodeAttributes(), net.getConditions());
		resultNet.setNodesAndEdges(net.getNodes(), net.getEdges());
		fullCleaning(resultNet, eo, task, false, true);

		return resultNet;
	}

	/**
	 * Clean an input reference network:
	 * Per pair of nodes, group all input edges into subclasses per root category of the EdgeOntology, unify the directionality
	 * (either all symmetric or all directed, as dicated by the edge ontology), and resolve conflicts within a root category.
	 * 
	 * @param net the network that needs cleaning
	 * @param eo the edge ontology
	 * @param task the task object that keeps track of the progress of this calculation (can be null)
	 * 
	 * @return a cleaned reference network representing the same semantic information
	 */
	public ReferenceNetwork fullInputRefCleaning(ReferenceNetwork net, EdgeOntology eo, ScheduledTask task)
	{
		ReferenceNetwork resultNet = new ReferenceNetwork(net.getName(), net.getID(), net.getAllNodeAttributes());
		resultNet.setNodesAndEdges(net.getNodes(), net.getEdges());
		fullCleaning(resultNet, eo, task, true, true);
		return resultNet;
	}

	/**
	 * Clean a generic input network:
	 * Per pair of nodes, group all input edges into subclasses per root category of the EdgeOntology, unify the directionality
	 * (either all symmetric or all directed, as dictated by the edge ontology), and resolve conflicts within a root category.
	 * 
	 * @param net the network that needs cleaning
	 * @param eo the edge ontology
	 * @param task the task object that keeps track of the progress of this calculation (can be null)
	 * 
	 * @return a cleaned reference network representing the same semantic information
	 */
	public InputNetwork fullInputCleaning(InputNetwork net, EdgeOntology eo, ScheduledTask task)
	{
		if (net instanceof ReferenceNetwork)
		{
			return fullInputRefCleaning((ReferenceNetwork) net, eo, task);
		}
		if (net instanceof ConditionNetwork)
		{
			return fullInputConditionCleaning((ConditionNetwork) net, eo, task);
		}
		
		InputNetwork resultNet = new InputNetwork(net.getName(), net.getID(), net.getAllNodeAttributes());
		resultNet.setNodesAndEdges(net.getNodes(), net.getEdges());
		fullCleaning(resultNet, eo, task, false, true);

		return resultNet;
	}

	/**
	 * For each node pair and each semantic root category, resolve conflicts by calling 
	 * {@link #resolveToOnePerRoot(Network, EdgeOntology, Set, Set, Node, Node)}.
	 * 
	 * This method is best run after {@link #removeRedundantEdges} as this already resolves some obvious conflicts.
	 * Further, it only works for input/consensus networks as it uses the source categories of the edge ontology.
	 * 
	 * @param net the network that needs cleaning
	 * @param eo the edge ontology
	 */
	protected void resolveToOne(Network net, EdgeOntology eo)
	{
		Set<Node> allNodes = net.getNodes();
		Set<Edge> newEdges = new HashSet<Edge>();
		Set<String> roots = eo.retrieveAllSourceRootCats(true);

		// first, determine all node pairs which are relevant in this network
		Set<Edge> oldEdges = net.getEdges();
		Map<String, Set<String>> directed_pairs = new HashMap<String, Set<String>>();
		Map<String, Set<String>> symmetrical_pairs = new HashMap<String, Set<String>>();

		Map<String, Node> mappedNodes = new HashMap<String, Node>();

		for (Edge e : oldEdges)
		{
			Node source = e.getSource();
			Node target = e.getTarget();

			String sourceID = source.getID();
			String targetID = target.getID();

			// edge is directed: store source -> target as a valid pair in the network
			if (!e.isSymmetrical())
			{
				if (!directed_pairs.containsKey(sourceID))
				{
					directed_pairs.put(sourceID, new HashSet<String>());
				}
				directed_pairs.get(sourceID).add(targetID);
			}

			// if the edge is symmetrical, store X -> Y with X's ID smaller than Y's ID
			if (e.isSymmetrical())
			{
				if (sourceID.compareTo(targetID) < 0)
				{
					if (!symmetrical_pairs.containsKey(sourceID))
					{
						symmetrical_pairs.put(sourceID, new HashSet<String>());
					}
					symmetrical_pairs.get(sourceID).add(targetID);
				}
				else
				{
					if (!symmetrical_pairs.containsKey(targetID))
					{
						symmetrical_pairs.put(targetID, new HashSet<String>());
					}
					symmetrical_pairs.get(targetID).add(sourceID);
				}
			}

			mappedNodes.put(sourceID, source);
			mappedNodes.put(targetID, target);
		}

		// For each node pair, perform a cleaning step for the asymmetrical pairs
		for (String ID1 : directed_pairs.keySet())
		{
			Node n1 = mappedNodes.get(ID1);
			Set<String> partners = directed_pairs.get(ID1);

			for (String ID2 : partners)
			{
				Node n2 = mappedNodes.get(ID2);

				Set<EdgeDefinition> edges = net.getAllEdgeDefinitions(n1, n2, false);
				if (!edges.isEmpty())
				{
					Map<String, Set<EdgeDefinition>> cleanedEdges = resolveToOnePerRoot(net, eo, roots, edges, n1, n2);
					for (Set<EdgeDefinition> cleans : cleanedEdges.values())
					{
						for (EdgeDefinition def : cleans)
						{
							newEdges.add(new Edge(n1, n2, def));
						}
					}
				}
			}
		}
		
		// For each node pair, perform a cleaning step for the symmetrical pairs
		for (String ID1 : symmetrical_pairs.keySet())
		{
			Node n1 = mappedNodes.get(ID1);
			Set<String> partners = symmetrical_pairs.get(ID1);

			for (String ID2 : partners)
			{
				Node n2 = mappedNodes.get(ID2);

				Set<EdgeDefinition> edges = net.getAllEdgeDefinitions(n1, n2, true);
				if (!edges.isEmpty())
				{
					Map<String, Set<EdgeDefinition>> cleanedEdges = resolveToOnePerRoot(net, eo, roots, edges, n1, n2);
					for (Set<EdgeDefinition> cleans : cleanedEdges.values())
					{
						for (EdgeDefinition def : cleans)
						{
							newEdges.add(new Edge(n1, n2, def));
						}
					}
				}
			}
		}
		net.setNodesAndEdges(allNodes, newEdges);
	}

	/**
	 * Clean a network:
	 * Group all input edges into subclasses per root category of the EdgeOntology, resolve conflicts within a root category.
	 * 
	 * This method only works for input/consensus networks as it uses the source categories of the edge ontology.
	 * 
	 * @param net the network that needs cleaning
	 * @param eo the edge ontology
	 * @param roots the root categories of the edge ontology
	 * @param oldEdges all edges between two nodes, including both directed and symmetrical edges
	 * @param source the source node (used for logging)
	 * @param target the target node (used for logging)
	 * 
	 * @return all edges grouped by semantic root category, with unified directionality, and only 1 edge per network and root cat.
	 */
	
    protected Map<String, Set<EdgeDefinition>> resolveToOnePerRoot(Network net, EdgeOntology eo, Set<String> roots, Set<EdgeDefinition> oldEdges, 
    		Node source, Node target)
	{
		Map<String, Set<EdgeDefinition>> mappedNormalEdges = getEdgesPerRoot(eo, roots, oldEdges);
		Map<String, Set<EdgeDefinition>> mappedSingleEdges = new HashMap<String, Set<EdgeDefinition>>();

		for (String rootCat : mappedNormalEdges.keySet())
		{
			Set<EdgeDefinition> normalEdges = mappedNormalEdges.get(rootCat);
			Set<EdgeDefinition> singleEdges = getMaxWeightedEdge(normalEdges, eo, net.getName(), source, target, rootCat);
			mappedSingleEdges.put(rootCat, singleEdges);
		}
		return mappedSingleEdges;
	}
	

	/**
	 * Group all input edges into subclasses per root category of the EdgeOntology.
	 * 
	 * @param eo the edge ontology
	 * @param roots the root categories of the edge ontology
	 * @param oldEdges the original set of input edges
	 * @return all input edges grouped by edge root category
	 */
	protected Map<String, Set<EdgeDefinition>> getEdgesPerRoot(EdgeOntology eo, Set<String> roots, Set<EdgeDefinition> oldEdges)
	{
		Map<String, Set<EdgeDefinition>> mappedEdges = new HashMap<String, Set<EdgeDefinition>>();

		for (EdgeDefinition refE : oldEdges)
		{
			String edgeType = refE.getType();

			int foundRoot = 0;

			for (String root : roots)
			{
				boolean belongsToRoot = eo.isSourceTypeChildOf(edgeType, root) >= 0;
				if (belongsToRoot)
				{
					foundRoot++;
					if (!mappedEdges.containsKey(root))
					{
						mappedEdges.put(root, new HashSet<EdgeDefinition>());
					}
					mappedEdges.get(root).add(refE);
				}
			}
			if (foundRoot == 0)
			{
				String errorMsg = " Edge source type '" + edgeType + "' could not be linked to a semantic root category in the edge ontology";
				throw new IllegalArgumentException(errorMsg);
			}
			if (foundRoot > 1)
			{
				String errorMsg = " Edge source type '" + edgeType + "' could be linked to more than one semantic root category in the edge ontology";
				throw new IllegalArgumentException(errorMsg);
			}
		}

		return mappedEdges;
	}

	/**
	 * Return the edge with the highest weight, or, if there are more, the most specific one.
	 * Affirmative and negated edges are not compared to eachother, and thus the highest one of both classes can be maintained.
	 * 
	 * It is assumed that getEdgesPerRoot was previously used to provide a set of edges which only contains edges for one root category,
	 * and that all edges within this category are either symmetrical, or all directed.
	 * 
	 * This method only works for input/consensus networks as it uses the source categories of the edge ontology.
	 * 
	 * @param edges the original set of input edges
	 * @param eo the edge ontology
	 * @param network_name the name/type of input network (used for logging)
	 * @param source the source node (used for logging)
	 * @param target the target node (used for logging)
	 * @param rootCat the root category (used for logging)
	 * 
	 * @return one edge, produced after resolving conflicts using the edge weights, or throws a RunTimeException if no best edge could be found.
	 */
	protected Set<EdgeDefinition> getMaxWeightedEdge(Set<EdgeDefinition> edges, EdgeOntology eo, String network_name, Node source, Node target, String rootCat)
	{
		double maxAffWeight = Double.NEGATIVE_INFINITY;
		double maxNegWeight = Double.NEGATIVE_INFINITY;
		
		boolean affSymm = true;
		boolean negSymm = true;
		
		int nr_aff = 0;
		int nr_neg = 0;

		for (EdgeDefinition e : edges)
		{
			double weight = e.getWeight();
			if (e.isNegated())
			{
				nr_neg++;
				maxNegWeight = Math.max(maxNegWeight, weight);
				if (! e.isSymmetrical())
				{
					negSymm = false;
				}
			}
			else
			{
				nr_aff++;
				maxAffWeight = Math.max(maxAffWeight, weight);
				if (! e.isSymmetrical())
				{
					affSymm = false;
				}
			}
		}

		Set<EdgeDefinition> thickestEdges = new HashSet<EdgeDefinition>();

		Set<String> affirmative_cats = new HashSet<String>();
		Set<String> negated_cats = new HashSet<String>();

		for (EdgeDefinition e : edges)
		{
			if (e.isNegated())
			{
				if (Math.abs(maxNegWeight - e.getWeight()) < 0.000001)		// allow for small rounding errors
				{
					thickestEdges.add(e);
					String type = eo.getSourceCategory(e.getType());
					negated_cats.add(type);
				}
			}
			else
			{
				if (Math.abs(maxAffWeight - e.getWeight()) < 0.000001)		// allow for small rounding errors
				{
					thickestEdges.add(e);
					String type = eo.getSourceCategory(e.getType());
					affirmative_cats.add(type);
				}
			}
		}

		// For the non-negated types, we take the most specific edge category that still covers all 
		String affirmativeConsensus = null;
		for (String aff_cat : affirmative_cats)
		{
			if (affirmativeConsensus == null)
			{
				affirmativeConsensus = aff_cat;
			}
			else
			{
				int child = eo.isSourceCatChildOf(aff_cat, affirmativeConsensus);   // if positive, aff_type is a child of the consensus
				int parent = eo.isSourceCatChildOf(affirmativeConsensus, aff_cat);  // if positive, aff_type is a parent of the consensus

				// they are siblings or something such: take the first common parent
				if (child < 0 && parent < 0)
				{
					Set<String> cats = new HashSet<String>();
					cats.add(aff_cat);
					cats.add(affirmativeConsensus);
					affirmativeConsensus = eo.retrieveFirstCommonParent(cats);
				}
				if (child > 0)
				{
					// this one is more specific!
					affirmativeConsensus = aff_cat;
				}
			}
		}

		// For the negated types, we take the most general one
		String negatedConsensus = null;
		for (String neg_cat : negated_cats)
		{
			if (negatedConsensus == null)
			{
				negatedConsensus = neg_cat;
			}
			else
			{
				int child = eo.isSourceCatChildOf(neg_cat, negatedConsensus);   // if positive, neg_type is a child of the consensus
				int parent = eo.isSourceCatChildOf(negatedConsensus, neg_cat);  // if positive, neg_type is a parent of the consensus

				// they are siblings or something such: take the first common parent
				// This is not entirely correct because negative evidence shouldn't travel up the tree... 
				// but it seems the most sensible thing to do to summarize the given information
				if (child < 0 && parent < 0)
				{
					Set<String> cats = new HashSet<String>();
					cats.add(neg_cat);
					cats.add(negatedConsensus);
					negatedConsensus = eo.retrieveFirstCommonParent(cats);
				}
				if (parent > 0)
				{
					// this one is more general!
					negatedConsensus = neg_cat;
				}
			}
		}

		Set<EdgeDefinition> results = new HashSet<EdgeDefinition>();
		
		if (affirmativeConsensus != null)
		{
			EdgeDefinition ed = new EdgeDefinition(affirmativeConsensus, affSymm, maxAffWeight, false);
			results.add(ed);
			
			if (nr_aff > 1)
			{
				logger.log("  Kept only the affirmative edge with weight (" + maxAffWeight + ") and type " + affirmativeConsensus + " between " + source + " and " + target + " for the category " + rootCat + " in " + network_name);
			}
		}
		if (affirmativeConsensus == null && nr_aff > 0)
		{
			String errorMsg = "Could not resolve the set of affirmative edges to one !";
			throw new IllegalArgumentException(errorMsg);
		}

		if (negatedConsensus != null)
		{
			EdgeDefinition ed = new EdgeDefinition(negatedConsensus, negSymm, maxNegWeight, true);
			results.add(ed);
			
			if (nr_neg > 1)
			{
				logger.log("  Kept only the negated edge with weight (" + maxNegWeight + ") and type " + negatedConsensus + " between " + source + " and " + target + " for the category " + rootCat + " in " + network_name);
			}
		}
		if (negatedConsensus == null && nr_neg > 0)
		{
			String errorMsg = "Could not resolve the set of negated edges to one !";
			throw new IllegalArgumentException(errorMsg);
		}

		return results;
	}

}
