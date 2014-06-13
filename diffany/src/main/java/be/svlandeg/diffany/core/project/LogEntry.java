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
	 * TODO
	 * @param date
	 * @param message
	 */
	public LogEntry(Date date, String message)
	{
		this.date = date;
		this.message = message;
	}
	
	/**
	 * TODO
	 * @return
	 */
	public String getMessage()
	{
		return message;
	}
	
	/**
	 * TODO
	 * @return
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
