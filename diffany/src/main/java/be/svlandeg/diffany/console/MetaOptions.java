package be.svlandeg.diffany.console;

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


import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

/**
 * This class defines the various 'meta' options available when running the Diffany jar through the console,
 * such as the help message and the version number.
 * 
 * @author Sofie Van Landeghem
 */
public class MetaOptions
{

	private Options options;
	
	protected static String helpShort = "h";
	protected static String versionShort = "v";

	/**
	 * Constructor initializes the meta options for this program
	 */
	public MetaOptions()
	{
		defineOptions();
	}
	
	/**
	 * Retrieve the options object for the meta data
	 * @return the options object
	 */
	public Options getMetaOptions()
	{
		return options;
	}

	
	/**
	 * Define the options available for providing meta data (all simple flags)
	 */
	private void defineOptions()
	{
		boolean hasArgument = false;
		
		options = new Options();
		
		options.addOption(new Option(versionShort, "version", hasArgument, "print the version information and exit"));
		options.addOption(new Option(helpShort, "help", hasArgument, "print this help message"));
	}
	
}
