package be.svlandeg.diffany.cytoscape;

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

import org.cytoscape.session.events.SessionAboutToBeSavedEvent;
import org.cytoscape.session.events.SessionAboutToBeSavedListener;
import org.cytoscape.session.events.SessionLoadedEvent;
import org.cytoscape.session.events.SessionLoadedListener;

/**
 * A listener that catches session load and save events from Cytoscape.
 * 
 * At this point, this only makes sure the Diffany visual styles are preserved
 * when the user loads a new session.
 * 
 * @author Thomas Van Parys
 *
 */
public class SessionListener implements SessionAboutToBeSavedListener,
		SessionLoadedListener {

	private Model model;
	
	
	/**
	 * Construct new listener
	 * @param model the Diffany {@link Model}
	 */
	public SessionListener(Model model) {
		this.model = model;
	}

	@Override
	public void handleEvent(SessionLoadedEvent arg0) {
		//re-register Visual Styles
		model.getSourceStyle().registerStyle();
		model.getDiffStyle().registerStyle();
	}

	@Override
	public void handleEvent(SessionAboutToBeSavedEvent arg0) {
		// TODO save CyProject information

	}

}
