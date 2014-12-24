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

import org.cytoscape.model.CyNetwork;
import org.cytoscape.work.TaskFactory;
import org.cytoscape.work.TaskIterator;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;

/**
 * Factory to run the {@link UpdateVisualStyleTask}
 * 
 * @author Thomas Van Parys
 *
 */
public class UpdateVisualStyleTaskFactory implements TaskFactory {

	private Model model;
	private CyProject cyProject;
	
	/**
	 * Construct a new factory to create {@link UpdateVisualStyleTask}s
	 * 
	 * @param model Diffany {@link Model}
	 * @param project the {@link CyProject} containing the used {@link CyNetwork}s.
	 */
	public UpdateVisualStyleTaskFactory(Model model, CyProject project) {
		this.model = model;
		this.cyProject = project;
	}

	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator it = new TaskIterator(new UpdateVisualStyleTask(model, cyProject));
		return it;
	}

	@Override
	public boolean isReady() {
		return true;
	}

}
