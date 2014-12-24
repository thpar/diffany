package be.svlandeg.diffany.core.networks;

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

import be.svlandeg.diffany.core.io.EdgeIO;

/**
 * Class that represents an edge in a network: an edge has a source and target node
 * and can have a certain weight. It is symmetrical or not, and may or may not be negated.
 * An edge can not be altered.
 * 
 * Edge data can be saved and loaded through the {@link EdgeIO} class
 * 
 * @author Sofie Van Landeghem
 */
public class Edge 
{
	protected Node source;
	protected Node target;
	protected EdgeDefinition def;

	/**
	 * Create a new node from a certain definition and specifying source and target nodes.
	 * The EdgeDefinition object is kept as such (not copied).
	 * 
	 * @param source the source node
	 * @param target the target node
	 * @param def the edge definition specifying the type, weight, symmetry and negation of the edge
	 */
	public Edge(Node source, Node target, EdgeDefinition def)
	{
		if (source == null)
		{
			String errormsg = "The source node should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		if (target == null)
		{
			String errormsg = "The target node should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.source = source;
		this.target = target;
		this.def = def;
	}

	/**
	 * Create a new edge with specified source and target nodes, direction and weight.
	 * The input parameters will be used to create an internal EdgeDefinition object.
	 * 
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param weight the weight or confidence of this edge (should be positive)
	 * @param negated defines whether or not the edge is negated (e.g. does NOT bind)
	 * @throws IllegalArgumentException when the specified weight is a negative number
	 */
	public Edge(String type, Node source, Node target, boolean symmetrical, double weight, boolean negated) throws IllegalArgumentException
	{
		this(source, target, new EdgeDefinition(type, symmetrical, weight, negated));
	}

	/**
	 * Create a new edge with default weight of 1.
	 * 
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param negated defines whether or not the edge is negated (e.g. does NOT bind)
	 */
	public Edge(String type, Node source, Node target, boolean symmetrical, boolean negated)
	{
		this(type, source, target, symmetrical, EdgeDefinition.DEFAULT_WEIGHT, negated);
	}

	/**
	 * Create a new edge, which will be defined as not negated.
	 * 
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 * @param weight the weight or confidence of this edge (should be positive)
	 */
	public Edge(String type, Node source, Node target, boolean symmetrical, double weight)
	{
		this(type, source, target, symmetrical, weight, EdgeDefinition.DEFAULT_NEG);
	}

	/**
	 * Create a new edge with default weight of 1.0 and negation off.
	 * 
	 * @param type the interaction type of this edge
	 * @param source the source node
	 * @param target the target node
	 * @param symmetrical defines whether the edge is symmetrical or directed from source to target
	 */
	public Edge(String type, Node source, Node target, boolean symmetrical)
	{
		this(type, source, target, symmetrical, EdgeDefinition.DEFAULT_WEIGHT, EdgeDefinition.DEFAULT_NEG);
	}

	/**
	 * Get the source node of this edge.
	 * When the edge is symmetrical, source and target will be defined consistently,
	 * though they can be used interchangeably in upstream code.
	 * 
	 * @return the source node of this edge
	 */
	public Node getSource()
	{
		return source;
	}

	/**
	 * Get the target node of this edge.
	 * When the edge is symmetrical, source and target will be defined consistently,
	 * though they can be used interchangeably in upstream code.
	 * 
	 * @return the target node of this edge
	 */
	public Node getTarget()
	{
		return target;
	}

	@Override
	public String toString()
	{
		return EdgeIO.writeToTab(this);
	}
	
	/**
	 * Get the definition of this edge
	 * @return the edge definition
	 */
	public EdgeDefinition getDefinition()
	{
		return def;
	}
	
	/**
	 * Get the type of this edge
	 * @return the edge type
	 */
	public String getType()
	{
		return def.type;
	}
	
	/**
	 * Return whether or not this edge is symmetrical 
	 * (if not, it goes specifically from source to target, otherwise there is no real direction)
	 * @return whether or not this edge is symmetrical
	 */
	public boolean isSymmetrical()
	{
		return def.symmetrical;
	}
	
	/**
	 * Return whether or not this edge is negated (e.g. does NOT bind) 
	 * @return whether or not this edge is negated
	 */
	public boolean isNegated()
	{
		return def.negated;
	}
	
	/**
	 * Get the weight of this edge.
	 * When it never has been set, it is 1.0 by default.
	 * The edge weight should always be positive. When it is 0, the edge can be ignored.
	 * @return the edge weight
	 */
	public double getWeight()
	{
		return def.weight;
	}

}
