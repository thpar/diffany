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
import org.cytoscape.work.TaskManager;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.tasks.UpdateVisualStyleTask;
import be.svlandeg.diffany.cytoscape.tasks.UpdateVisualStyleTaskFactory;

/**
 * Run the {@link UpdateVisualStyleTask}
 * 
 * @author Thomas Van Parys
 *
 */
public class UpdateVisualStyleAction extends AbstractCyAction {
	
	private static final long serialVersionUID = -8382501915601785381L;
	private static final String BUTTON_TITLE = "Update Visual Style";
	private Model model;
	private CyProject cyProject;
	
	/**
	 * Action that refreshes the registered visual styles and reapplies them on the
	 * Views, executing the {@link UpdateVisualStyleTask}
	 * 
	 * @param model the Diffany {@link Model}
	 */
	public UpdateVisualStyleAction(Model model){
		super(BUTTON_TITLE);
		this.model = model;
	}
	
	/**
	 * 
	 * @param model the Diffany {@link Model}
	 * @param cyProject {@link CyProject} containing the {@link CyNetwork}s to update.
	 */
	public UpdateVisualStyleAction(Model model, CyProject cyProject){
		super(BUTTON_TITLE);
		this.model = model;
		this.cyProject = cyProject;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		UpdateVisualStyleTaskFactory tf;
		if (cyProject != null){
			tf = new UpdateVisualStyleTaskFactory(model, cyProject);			
		} else {
			tf = new UpdateVisualStyleTaskFactory(model, model.getSelectedProject());			
		}

		if (tf.isReady()){
			TaskIterator it = tf.createTaskIterator();			
			TaskManager<?, ?> tm = model.getServices().getTaskManager();
			tm.execute(it);
		}

	}

}
