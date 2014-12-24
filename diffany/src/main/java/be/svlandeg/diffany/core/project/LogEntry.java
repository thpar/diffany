package be.svlandeg.diffany.core.project;

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
 * A LogEntry is a logged message together with relevant meta data, used by the Logger object.
 * 
 * @author Sofie Van Landeghem
 */
public class LogEntry
{
	
	private Date date;
	private String message;
	
	/**
	 * Create a new log entry at a specific date and with a specific message
	 * @param date the date at which the log entry occurred
	 * @param message a specification of the logged events
	 */
	public LogEntry(Date date, String message)
	{
		this.date = date;
		this.message = message;
	}
	
	/**
	 * Retrieve the message of this log entry
	 * @return the message of this log entry
	 */
	public String getMessage()
	{
		return message;
	}
	
	/**
	 * Retrieve the date of this log entry
	 * @return the date of this log entry
	 */
	public Date date()
	{
		return date;
	}
	
	public String toString()
	{
		return date + " : " + message;
	}

}
