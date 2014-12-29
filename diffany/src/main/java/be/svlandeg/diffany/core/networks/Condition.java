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

import java.util.HashSet;
import java.util.Set;

/**
 * This class describes an experimental condition that can be linked to a certain network.
 * A condition is described in free text, but may also be associated to ontology terms.
 * 
 * @author Sofie Van Landeghem
 *
 */
public class Condition
{

	protected String description;
	protected Set<String> ontologies;

	/**
	 * Create a new condition that can be linked to a network.
	 * 
	 * @param description free-text description of the condition
	 * @param ontologies set of corresponding ontology terms (should not be null)
	 * 
	 * @throws IllegalArgumentException when the description or the set of ontologies is null
	 */
	public Condition(String description, Set<String> ontologies) throws IllegalArgumentException
	{
		if (description == null)
		{
			String errormsg = "The description should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.description = description;
		setOntologyTerms(ontologies);
	}

	/**
	 * Create a new condition that can be linked to a network.
	 * The set of ontology terms will be initialized to an empty set.
	 * 
	 * @param description free-text description of the condition
	 * @throws IllegalArgumentException when the description is null
	 */
	public Condition(String description) throws IllegalArgumentException
	{
		this(description, new HashSet<String>());
	}

	/**
	 * Cloning constructor
	 * @param oldCondition the condition that needs cloning
	 */
	public Condition(Condition oldCondition)
	{
		this(new String(oldCondition.description), new HashSet<String>(oldCondition.ontologies));
	}

	/**
	 * Return a free-text description of the condition
	 * @return the description of the condition
	 */
	public String getDescription()
	{
		return description;
	}

	/**
	 * Return the corresponding ontology terms for this condition
	 * @return the ontology terms for this condition (may be empty but should not be null)
	 */
	public Set<String> getOntologyTerms()
	{
		return ontologies;
	}

	/**
	 * Define the set of ontologies corresponding to this condition 
	 * (overwrites all previously defined terms).
	 * 
	 * @param ontologies the set of ontology terms (should not be null)
	 * @throws IllegalArgumentException when the description or the set of ontologies is null
	 */
	public void setOntologyTerms(Set<String> ontologies) throws IllegalArgumentException
	{
		if (ontologies == null)
		{
			String errormsg = "The set of ontologies should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.ontologies = ontologies;
	}

	/**
	 * Add a number of ontology terms to the set of previously defined terms.
	 * @param ontology a new ontology term
	 */
	public void addOntologyTerm(String ontology)
	{
		ontologies.add(ontology);
	}

	public String toString()
	{
		return description;
	}
}
