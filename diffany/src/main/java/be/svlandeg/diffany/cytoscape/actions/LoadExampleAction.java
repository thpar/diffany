package be.svlandeg.diffany.cytoscape.actions;

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

import java.awt.event.ActionEvent;

import org.cytoscape.application.swing.AbstractCyAction;
import org.cytoscape.model.CyNetwork;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.swing.DialogTaskManager;

import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.cytoscape.internal.Services;
import be.svlandeg.diffany.cytoscape.tasks.LoadExampleTaskFactory;
import be.svlandeg.diffany.examples.GenericExample;

/**
 * Action that loads example {@link CyNetwork}s into a Cytoscape session, based on a given {@link Project}.
 * 
 * @author Thomas Van Parys
 *
 */
public class LoadExampleAction extends AbstractCyAction{

	private static final long serialVersionUID = 1L;
	private Services services;
	private Project exampleProject;
	private int runConfigurationID;
	private GenericExample example;
	
	/**
	 * Load an example {@link Project} into the current Cytoscape session. 
	 * Project settings are ignored. Only the input networks are read and translated into {@link CyNetwork}s.
	 * 
	 * @param services the Cytoscape {@link Services}s
	 * @param name display name in the menu
	 * @param exampleProject {@link Project} to load the source networks from
	 * @param runConfigurationID id of the example configuration to run in the project
	 */
	public LoadExampleAction(Services services, String name, Project exampleProject, int runConfigurationID) {
		super(name, services.getCyApplicationManager(), null, null);
		this.services = services;
		this.exampleProject = exampleProject;
		this.runConfigurationID = runConfigurationID;
		setPreferredMenu("Apps.Diffany.Examples");
	}

	/**
	 * Set up the action, but don't load the example just yet. The project will only be constructed upon
	 * runtime.
	 * 
	 * @param services the services
	 * @param example the example
	 */
	public LoadExampleAction(Services services, GenericExample example) {
		super(example.getName(), services.getCyApplicationManager(), null, null);
		this.services = services;
		this.example = example;
		setPreferredMenu("Apps.Diffany.Examples");
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		LoadExampleTaskFactory tf = null;
		if (exampleProject == null){
			tf = new LoadExampleTaskFactory(services, example);			
		} else {
			tf = new LoadExampleTaskFactory(services, exampleProject, runConfigurationID);			
		}
		

		if (tf.isReady()){
			TaskIterator it = tf.createTaskIterator();			
			DialogTaskManager dtm = services.getDialogTaskManager();
			dtm.execute(it);
		}

	}

}
