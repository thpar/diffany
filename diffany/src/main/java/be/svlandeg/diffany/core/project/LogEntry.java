package be.svlandeg.diffany.core.project;

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
