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

import be.svlandeg.diffany.cytoscape.Model;

/**
 * Factory to run {@link RunProjectTask}, followed by a {@link UpdateVisualStyleTask}
 * 
 * @author Thomas Van Parys
 *
 */
public class RunProjectTaskFactory implements TaskFactory {

	private Model model;
	
	/**
	 * Construct a new factory.
	 * 
	 * @param model Diffany {@link Model}
	 */
	public RunProjectTaskFactory(Model model) {
		this.model = model;
	}

	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator it = new TaskIterator();
		
		it.append(new RunProjectTask(model, model.getSelectedProject()));
		it.append(new UpdateVisualStyleTask(model, model.getSelectedProject()));
		
		return it;
	}

	@Override
	public boolean isReady() {
		return model.getSelectedProject().canExecute(model);
	}

}
