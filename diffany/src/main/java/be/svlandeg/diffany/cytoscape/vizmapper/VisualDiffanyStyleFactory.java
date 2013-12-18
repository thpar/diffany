package be.svlandeg.diffany.cytoscape.vizmapper;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.internal.Services;

public class VisualDiffanyStyleFactory {
	
	public enum Type{
		SOURCE, DIFF;
	}
	
	static public void registerNewVisualStyle(Type type, Model model){
		Services services = model.getServices();
		switch(type){
		case SOURCE:	
			VisualSourceStyle sourceStyle = new VisualSourceStyle(services);
			model.getGuiModel().setSourceStyle(sourceStyle);
			break;
		case DIFF:
			VisualDiffStyle diffStyle = new VisualDiffStyle(services);
			model.getGuiModel().setDiffStyle(diffStyle);
			break;
		}
	}
}
