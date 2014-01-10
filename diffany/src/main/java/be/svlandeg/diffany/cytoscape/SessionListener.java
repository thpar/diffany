package be.svlandeg.diffany.cytoscape;

import org.cytoscape.session.events.SessionAboutToBeSavedEvent;
import org.cytoscape.session.events.SessionAboutToBeSavedListener;
import org.cytoscape.session.events.SessionLoadedEvent;
import org.cytoscape.session.events.SessionLoadedListener;

import be.svlandeg.diffany.cytoscape.vizmapper.VisualDiffStyle;

public class SessionListener implements SessionAboutToBeSavedListener,
		SessionLoadedListener {

	private Model model;
	
	
	
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
