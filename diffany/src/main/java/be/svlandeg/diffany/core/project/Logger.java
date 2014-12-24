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

import java.util.ArrayList;
import java.util.Date;
import java.util.List;

/**
 * A Logger object records all messages that are relevant to the user of the networks algorithms.
 * For instance, it reports on how edge conflicts were resolved, when default values were used, etc.
 * 
 * @author Sofie Van Landeghem
 */
public class Logger
{
	private List<LogEntry> logs;

	/**
	 * Create a new instance with an empty log file.
	 */
	public Logger()
	{
		clean();
	}

	/**
	 * Clean all previous logs.
	 */
	public void clean()
	{
		logs = new ArrayList<LogEntry>();
	}

	/**
	 * Retrieve all the logged messages.
	 * 
	 * @return all logged messages
	 */
	public List<LogEntry> getAllLogMessages()
	{
		return logs;
	}

	/**
	 * Record a log message, with the date specified to be at the time of this function call.
	 * 
	 * @param message the log message
	 */
	public void log(String message)
	{
		logs.add(new LogEntry(new Date(), message));
	}

}
