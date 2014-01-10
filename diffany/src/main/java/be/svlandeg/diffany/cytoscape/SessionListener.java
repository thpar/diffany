package be.svlandeg.diffany.cytoscape;

import org.cytoscape.session.events.SessionAboutToBeSavedEvent;
import org.cytoscape.session.events.SessionAboutToBeSavedListener;
import org.cytoscape.session.events.SessionLoadedEvent;
import org.cytoscape.session.events.SessionLoadedListener;

/**
 * A listener that catches session load and save events from Cytoscape.
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
