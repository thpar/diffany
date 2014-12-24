package be.svlandeg.diffany.core.networks;

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
