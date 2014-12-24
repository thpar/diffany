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
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.swing.DialogTaskManager;

import be.svlandeg.diffany.cytoscape.CyProject;
import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.tasks.RunProjectTaskFactory;

/**
 * Action that runs the current {@link CyProject}.
 * 
 * @author Thomas Van Parys
 *
 */
public class RunProjectAction extends AbstractCyAction {

	private static final long serialVersionUID = 1L;
	private Model model;

	private static final String BUTTON_TITLE = "Start";
	
	/**
	 * Run the current {@link CyProject}
	 * 
	 * @param model the Diffany {@link Model}
	 * @param menuTitle title to use in the Apps menu
	 */
	public RunProjectAction(Model model, String menuTitle) {
		super(menuTitle, model.getServices().getCyApplicationManager(), null, null);
		setPreferredMenu("Apps.Diffany");
		this.model = model;
	}
	
	/**
	 * Run the current {@link CyProject}. As label, a default button label will be used.
	 * 
	 * @param model the Diffany {@link Model}
	 */
	public RunProjectAction(Model model){
		super(BUTTON_TITLE);
		this.model = model;
	}



	@Override
	public void actionPerformed(ActionEvent e) {
		RunProjectTaskFactory tf = new RunProjectTaskFactory(model);
		
		if (tf.isReady()){
			TaskIterator it = tf.createTaskIterator();			
			DialogTaskManager dtm = model.getServices().getDialogTaskManager();
			dtm.execute(it);
		}
	}

}
