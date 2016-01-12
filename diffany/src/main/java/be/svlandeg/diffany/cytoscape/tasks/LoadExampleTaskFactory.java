package be.svlandeg.diffany.cytoscape.tasks;

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

import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskIterator;

import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.cytoscape.internal.Services;
import be.svlandeg.diffany.examples.GenericExample;

/**
 * Factory to call {@link LoadExampleTask}
 * 
 * @author Thomas Van Parys
 *
 */
public class LoadExampleTaskFactory implements TaskFactory{

	private Services services;
	private Project exampleProject;
	private int runConfigurationID;
	private GenericExample example;
	
	/**
	 * 
	 * @param services the app {@link Services}
	 * @param exampleProject {@link Project} to be used as example input.
	 * @param runConfigurationID id of the example configuration to run in the project
	 */
	public LoadExampleTaskFactory(Services services, Project exampleProject, int runConfigurationID) {
		this.services = services;
		this.exampleProject = exampleProject;
		this.runConfigurationID = runConfigurationID;
	}
	
	/**
	 * 
	 * @param services the app {@link Services}
	 * @param example {@link Project} to be used as example input.
	 */
	public LoadExampleTaskFactory(Services services, GenericExample example) {
		this.services = services;
		this.example = example;
	}
	
	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator it = new TaskIterator();
		
		if (exampleProject != null){
			it.append(new LoadExampleTask(services, exampleProject, runConfigurationID));			
		} else {
			it.append(new LoadExampleTask(services, example));			
		}

		return it;
	}

	@Override
	public boolean isReady() {
		return true;
	}

}
