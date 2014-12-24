package be.svlandeg.diffany.core.progress;

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


import java.util.Date;


/**
 * This class implements the {@link ProgressListener} interface and prints the progress to 
 * the standard output stream, together with a time stamp.
 * 
 * @author Sofie Van Landeghem
 */
public class StandardProgressListener extends ProgressListener
{
	
	private boolean print;
	
	/**
	 * Constructor which specifies whether or not to print the progress messages.
	 * Usually, this value will be true, except for the JUnit usage of this class.
	 * @param print whether or not to print the progress messages
	 */
	public StandardProgressListener(boolean print)
	{
		super();
		this.print = print;
	}
	
	@Override
	protected void setProgress(String message, int progress, int total)
	{
		if (print)
		{
			System.out.println(new Date() + ": " + message + ": processed " + progress + " of " + total);
		}
	}
	
}
