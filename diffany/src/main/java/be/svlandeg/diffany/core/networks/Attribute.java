package be.svlandeg.diffany.core.networks;

public class Attribute {
	private String name;
	private Class<?> type;
	
	
	public Attribute(String name, Class<?> type){
		this.name = name;
		this.type = type;
	}


	@Override
	public boolean equals(Object obj) {
		if (obj instanceof Attribute){
			Attribute att = (Attribute)obj;
			return name.equals(att.name);
		}
		return false;
	}


	@Override
	public String toString() {
		return name;
	}


	public String getName() {
		return name;
	}


	public Class<?> getType() {
		return type;
	}
	
}
