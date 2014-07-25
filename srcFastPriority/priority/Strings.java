package priority;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * The class reads from a file called "strings"
 * all the strings that will be displayed as
 * part of the GUI.
 * @author raluca
 */
public class Strings {
	static public HashMap<String,String> strings = new HashMap<String,String>();
	static public String file_name = "./config/strings";
	
	/** Reads the strings from the configuration file and returns 
	 * an error message or the empty string. */
	static public String setStrings() {
		String[] result;
		String line;
		BufferedReader br = null;
		int line_no = 0;

		try {
			br = new BufferedReader(new FileReader(file_name));
			line = br.readLine();
			line_no++;
			
			while (line != null) {
				if (line.compareTo("") != 0) {
				     result = line.split("=");
				     if (result.length != 2) 
				    	return "Error: the configuration file " + file_name + " is not well formatted!\n" +
				    		 "(line " + line_no + ": \"" + line + "\")";
				     if (result[0].indexOf(' ') >= 0)
					    	return "Error: the configuration file " + file_name + " is not well formatted!\n" +
				    		 "(line " + line_no + ": \"" + line + "\")";				    	 
				     strings.put(result[0].trim(), result[1].trim());
				}
				line = br.readLine();
			}			
			
			br.close();
		}
		catch (IOException e) {
			try { br.close(); } 
			catch (Exception ee) { }
			return e.getMessage();
		}	
		return "";
	}
	
	/** Returns the string for the object objectName
	 * or "" if no such object exists. */
	public static String getString(String objectName) {
		if (strings.containsKey(objectName))
			return (String)strings.get(objectName);
		else return "";
	}
}
