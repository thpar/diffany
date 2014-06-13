package be.svlandeg.diffany.core.project;

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
	//TODO v2.0: use the Logger to contain Exceptions inheriting from SolvableException.


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
