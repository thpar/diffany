package be.svlandeg.diffany.cytoscape.vizmapper;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.internal.Services;

public class VisualDiffanyStyleFactory {
	
	public enum Type{
		SOURCE, DIFF;
	}
	
	static public void registerNewVisualStyle(Type type, Model model){
		VisualDiffanyStyle style = null;
		Services services = model.getServices();
		switch(type){
		case SOURCE:	
			style = new VisualSourceStyle(services);
			model.getGuiModel().setSourceStyle(style.getVisualStyle());
			break;
		case DIFF:
			style = new VisualDiffStyle(services);
			model.getGuiModel().setDiffStyle(style.getVisualStyle());
			break;
		}
	}
}
