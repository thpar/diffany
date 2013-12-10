package be.svlandeg.diffany.cytoscape.vizmapper;

import be.svlandeg.diffany.internal.Services;

public class VisualDiffanyStyleFactory {
	
	public enum Type{
		SOURCE, DIFF;
	}
	
	static public VisualDiffanyStyle registerNewVisualStyle(Type type, Services services){
		VisualDiffanyStyle style = null;
		switch(type){
		case SOURCE:	
			style = new VisualSourceStyle(services);
			break;
		case DIFF:
			style = new VisualDiffStyle(services);
			break;
		}
		return style;
	}
}
