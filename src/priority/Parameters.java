package priority;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.File;
import java.io.Serializable;
import java.util.HashMap;
import java.util.ArrayList;
import priority.gibbs.GibbsStatic;

/**
 * The class Parameters contains all the  parameters as static members that are available from all the other classes.
 * @author raluca
 * Some updates of code have been made by Carlos A. Sierra (carlos.andres.sierra.v@gmail.com)
 */
public class Parameters implements Serializable {
	private static final long serialVersionUID = 1;

	/* The hash that maps the var names to their initial values */
	static private HashMap<String, String> vars = new HashMap<String, String>();
	
	/* The file containing the initial values */
	static private String file_name = "./config/params"; 
	
	/* Path parameters: */	
	public static String fname_path;      /* the path for the FASTA files */
	public static String back_file;       /* the background file name */
	public static String path_output;     /* the path for the output files */
		
	/* Running parameters: */
	public static boolean noocflag;      /* "allow for no occurrence of the motif"-flag */
	public static boolean revflag;       /* flag for searching in the reverse strand also */
	public static boolean otherclass;    /* whether to allow for "other" class/ uniform prior */
	public static int trials;            /* the number of trials */
	public static int iter;              /* the number of iterations */
	public static int outputStep;        /* the output is printed after every outputStep iterations */
	public static boolean createLogl;    /* true=create the logl files */
	
	public static int wsize;             /* the window size */
	public static int wsizeMin;          /* the minimum window size */
	public static int wsizeMax;          /* the maximum window size */
	public static int wsizeMinFromPrior; /* the maximum window size (specified in a prior file) */
	public static int wsizeMaxFromPrior; /* the maximum window size (specified in a prior file) */

	public static int bkgrOrder;         /* the background model order */
	public static int bkgrOrderMin;      /* the minimum background model order */
	public static int bkgrOrderMax;      /* the maximum background model order */
	
	public static double phi_prior[] = new double[4];  /* dirichlet prior counts for the model parameters (phi) */
	public static double sampling_exp[][];             /* the exponents for prior and likeliood, used for sampling */ 
	
	/* Class/prior and TF variables */
	public static String prior_dirs[];  /* the class/prior directories (one dir for each class) */
	public static int putative_class;   /* the number of the class/prior the TF is likely to belong to */ 
	
	public static String tf_path;       /* the path to a file containing TF names */
	public static String tf_name;       /* a TF name (used when individualTF=true) */
	public static boolean individualTF; /* true = we use a single TF, false -> file with TF names */ 
	
		
	/* Parameters that are set exactly before starting the algorithm */
	/* (they are not read from the config file "params") */
	public static String tf_names[];
	public static double back[];
	
	/* Parameter used to check the background model */ 
	static private double precision;
	
	/* More running parameters */
	public static double d;
	public static double pseudocounts_class, pseudocounts_putative_class;
	public static double flat_prior_other_class;
	
	public static boolean multiple_priors; /* true = multiple prios (e.g.class), false = single/uniform prior */    
	
	/* Output parameters */
	public static int bestZ[] = null;
	public static int bestC[] = null;
	public static double bestPhi[][] = null;
	public static String comboseq[] = null;
	public static String comboseq_names[] = null;
	public static String outputTF_name = null;
	public static boolean outputParams_are_set = false;
	
	
	/**
	 * Reads the variables from the configuration file and returns an error message or the empty string.
	 * @return
	 */
	static public String setStrings() 
	{
		return Parameters.setStrings(Parameters.file_name);
	}

	
	/**
	 * Reads the variables from the given configuration file and returns an error message or the empty string.
	 * @param file_name
	 * @return error message or the empty string
	 */
	@SuppressWarnings("resource")
	static public String setStrings(String file_name) 
	{
		Parameters.file_name = file_name;
		String[] result;
		String line;
		BufferedReader br = null;
		int line_no = 0;

		try 
		{
			br = new BufferedReader(new FileReader(file_name));
			line = br.readLine();
			line_no++;
			
			while (line != null) 
			{
				if (line.compareTo("") != 0) 
				{
				     result = line.split("=");
				     
				     if (result.length != 2)
				     	return "Error: the configuration file " + file_name + " is not well formatted!\n" +
				    		 "(line " + line_no + ": \"" + line + "\")";
				     
				     if (result[0].indexOf(' ') >= 0)
				     	 return "Error: the configuration file " + file_name + " is not well formatted!\n" +
				    		 "(line " + line_no + ": \"" + line + "\")";
				     
				     vars.put(result[0].trim(), result[1].trim());
				}
				line = br.readLine();
			}			
			br.close();
		}
		catch (IOException e) 
		{
			try 
			{ 
				br.close(); 
			} 
			catch (Exception ee) { }
			
			return e.getMessage();
		}	
		return "";
	}

	
	/**
	 * Returns the string for the variable varName or "" if no such variable exists in the hash.
	 * @param varName
	 */
	public static String getValue(String varName) 
	{
		return (vars.containsKey(varName)) ? (String)vars.get(varName) : "";
	}
			 
	
	
	/**
	 * Sets the parameters to their default values, as given in the configuration file. 
	 * @param filename
	 * @return "" if successful and an error message if a value is not valid
	 */
	public static String setDefaultParameters(String filename) 
	{	
		Parameters.file_name = filename;
		return setDefaultParameters();
	}
	
	
	/**
	 * Sets the parameters to their default values, as given in the configuration file. 
	 * @return "" if successful and an error message if a value is not valid.
	 */
	public static String setDefaultParameters() 
	{	
		String err = setStrings();
		if(err.compareTo("") != 0)
			return err;
		
		/* otherclass/uniformprior it must be true or false */
		otherclass = ((getValue("otherpriortype")).compareTo("true") == 0) ? true : false;
		

	    /* it must have length >= 1 (otherwise it's just "otherclass/uniform" )
	     * and each entry must be the path of a readable directory */
	    prior_dirs = (getValue("prior_dirs")).split(";");
	    
		if(getValue("prior_dirs").compareTo("") == 0 ||	prior_dirs.length < 1) 
		{
			prior_dirs = new String[0];
	    	otherclass = true;
	    	multiple_priors = false;
		}
		else 
		{
			for(int i = 0; i < prior_dirs.length; i++) 
			{
				File dir = new File(prior_dirs[i]);
				
				if ((!dir.exists()) || (!dir.isDirectory()) || (!dir.canRead()))
					return "Error: the path for the prior files (\"" + prior_dirs[i] + "\") does not exist or it is not a readable directory!" +
		    			"\nPlease see the configuration file \"" + file_name + "\"!";
			}
		
			multiple_priors = ((prior_dirs.length == 1) && (otherclass == false)) ? false : true;
		}
		
				
		/* it must be between -1 and #classes-1 */
		String str = getValue("putative_priortype");
		try 
		{ 
			putative_class = Integer.parseInt(str); 
		}
		catch (Exception e) 
		{ 
			putative_class  = -1; 
		}
		
		if ((putative_class < -1) || (putative_class >= prior_dirs.length)) 
			return "Error: the value for the putative prior-type variable (" + putative_class + ") is not valid!\nIt must be between -1 and " + 
	    		(prior_dirs.length-1) + ".\nPlease see the configuration file \"" + file_name + "\"!";						
		
		
		/* *********************************************************************************** */
		/* it must be between 3 and 20 */
		str = getValue("wsizeMin");
		try 
		{ 
			wsizeMin = Integer.parseInt(str); 
		}
		catch (Exception e) 
		{ 
			wsizeMin = 3; 
		}
		
		
		if ((wsizeMin < 3) || (wsizeMin > 20)) 
			return "Error: the value for the wsizeMin variable (" +	wsizeMin + ") is not valid!\nIt must be between 3 and 20" + 
	    		".\nPlease see the configuration file \"" + file_name + "\"!";
		
		/* it must be between 3 and 20 and smaller than wsizeMin */
		str = getValue("wsizeMax");
		try 
		{ 
			wsizeMax = Integer.parseInt(str); 
		}
		catch (Exception e) 
		{ 
			wsizeMax = 20; 
		}
		
		
		if ((wsizeMax < wsizeMin) || (wsizeMax > 20))
			return "Error: the value for the wsizeMax variable (" + wsizeMax + ") is not valid!\nIt must be between " + wsizeMin + " and 20" + 
	    		".\nPlease see the configuration file \"" + file_name + "\"!";
		
		/* it must be between wsizeMin and wsizeMax */
		str = getValue("wsize");
		try 
		{ 
			wsize = Integer.parseInt(str); 
		}
		catch (Exception e) 
		{ 
			wsize = wsizeMin; 
		}
		
		
		if ((wsize < wsizeMin) || (wsize > wsizeMax))
	    	return "Error: the value for the wsize variable (" + wsizeMax + ") is not valid!\nIt must be between " + wsizeMin 
	    			+ " and " + wsizeMax +	".\nPlease see the configuration file \"" + file_name + "\"!";						
	
		/* *********************************************************************************** */
		/* they must be true or false */
		revflag = ((getValue("revflag")).compareTo("false") == 0) ? false : true;
		noocflag = ((getValue("noocflag")).compareTo("false") == 0) ? false : true;
		createLogl = ((getValue("logl")).compareTo("false") == 0) ? false : true;
		
		
		/* *********************************************************************************** */
		/* they must be >= 1 */
		str = getValue("trials");
		try 
		{ 
			trials = Integer.parseInt(str); 
		}
		catch (Exception e) 
		{ 
			trials = 1; 
		}
		
		if (trials < 1) 
			trials = 1;
				
		
		str = getValue("iter");
		try 
		{ 
			iter = Integer.parseInt(str); 
		}
		catch (Exception e)
		{ 
			iter = 1; 
		}
		
		if (iter < 1) 
			iter = 1;
				
		
		str = getValue("outputStep");
		try 
		{ 
			outputStep = Integer.parseInt(str); 
		}
		catch (Exception e) 
		{ 
			outputStep = 1; 
		}
		
		if (outputStep < 1) 
			outputStep = 1;
		
		
		str = getValue("precision");
		try 
		{ 
			precision = Double.parseDouble(str); 
		}
		catch (Exception e) 
		{ 
			precision = 0.0001; 
		}
		
		if ((precision < 0) || (precision >= 1)) 
			precision = 0.0001;
		

		/* *********************************************************************************** */
		/* it must be true or false */
		individualTF = ((getValue("individualTF")).compareTo("true") == 0) ? true : false;
		
		/* it must be a valid file or it can be "" if individualTF is true */
		tf_path = getValue("tf_path");
		File file = new File(tf_path);
		
		if ((!file.exists()) || (!file.isFile()) || (!file.canRead())) 
	    	if (individualTF == true)
	    		tf_path = "";
	    	else
	    		return "Error: the file containing TF names (\"" + tf_path + "\"),\nas specified " + "in the configuration file \"" + file_name + 
	    	       "\",\ndoes not exist or is not a readable file!";
		
		/* it must be a non-emtpy string if individualTF is true */
		tf_name = getValue("tf_name");
		if ((tf_name.compareTo("") == 0) && (individualTF == true))
			return "Error: the value for the TF name is not specified in the configuration file \"" + file_name + "\"!";
		
		/* *********************************************************************************** */
		/* they must be  valid directories */	    
		fname_path = getValue("fname_path"); 
		File dir = new File(fname_path);
	    if ((!dir.exists()) || (!dir.isDirectory()) || (!dir.canRead()))
	    	return "Error: the path for the FASTA files (\"" + fname_path + "\"),\nas specified in the configuration file \"" + file_name + "\",\n" + 
	    	       "does not exists or it is not a readable directory!";
	    	 
		path_output = getValue("path_output"); 
		dir = new File(path_output);
	    if ((!dir.exists()) || (!dir.isDirectory()) || (!dir.canWrite()))
	    	return "Error: the path for the output files (\"" + path_output + "\"),\nas specified in the configuration file \"" + file_name + "\",\n" +
	    	       "does not exists or it is not a writable directory!";
	    

	    /* *********************************************************************************** */
		/* it must be a valid file */
		back_file = getValue("back_file");
		file = new File(back_file);
		if ((!file.exists()) || (!file.isFile()) || (!file.canRead())) 
			return "Error: the background model file (\"" + back_file + "\"),\nas specified in the configuration file \"" + file_name + "\",\n" +
	    	       "does not exist or it is not a readable file!";		
		
		/* it must be >= 0 */
		str = getValue("bkgrOrderMin");
		try 
		{ 
			bkgrOrderMin = Integer.parseInt(str); 
		}
		catch (Exception e) 
		{ 
			bkgrOrderMin = 0; 
		}
		
		if (bkgrOrderMin < 0) 
			return "Error: the value for the bkgrOrderMin variable (" +	bkgrOrderMin + ") is not valid!\nIt must be > 0." + 
	    		" Please see the configuration file \"" + file_name + "\"!";						
			
		
		/* it must >= bkgrOrderMin */
		str = getValue("bkgrOrderMax");
		try 
		{ 
			bkgrOrderMax = Integer.parseInt(str); 
		}
		catch (Exception e) 
		{ 
			bkgrOrderMax = bkgrOrderMin; 
		}
		
		if (bkgrOrderMax < bkgrOrderMin)
			return "Error: the value for the bkgrOrderMax variable (" +	bkgrOrderMax + ") is not valid!\nIt must be > " + bkgrOrderMin + 
	    		". Please see the configuration file \"" + file_name + "\"!";						
		
		/* it must be between bkgrOrderMin and bkgrOrderMax */
		str = getValue("bkgrOrder");
		
		
		try 
		{ 
			bkgrOrder = Integer.parseInt(str); 
		}
		catch (Exception e) 
		{ 
			bkgrOrder = bkgrOrderMax; 
		}
		
		if ((bkgrOrder < bkgrOrderMin) || (bkgrOrder > bkgrOrderMax)) 
			return "Error: the value for the bkgrOrder variable (" + bkgrOrder + ") is not valid!\nIt must be between " + bkgrOrderMin + " and " + bkgrOrderMax + 
	    		".\nPlease see the configuration file \"" + file_name + "\"!\n";						
		
		/* *********************************************************************************** */ 
	    /* it must have length=4 and each entry must be a double positive value */
		String strs[] = (getValue("phi_prior")).split(" "); 
		if (strs.length != 4) 
			return "Error: the phi_prior variable specified in the configuration file \"" + file_name + "\" must contain exactly 4 real positive values!";
		
		for(int i = 0; i < 4; i++) 
		{
			try 
			{ 
				phi_prior[i] = Double.parseDouble(strs[i]); 
			}
			catch (Exception e) 
			{
		    	return "Error: one of the values in the phi_prior variable (" +  strs[i] +	"),\nas specified in the " +
		    		" configuration file \"" + file_name + "\", is not valid!\nThe phi_prior variable must contain only real positive values!";			
			}
			
			if (phi_prior[i] <= 0) 
				return "Error: one of the values in the phi_prior variable (" +  strs[i] + "),\nas specified in the " +
	    		" configuration file \"" + file_name + "\", is not valid!\nThe phi_prior variable must contain only real positive values!";
		}

		
		/* *********************************************************************************** */ 
	    /* it must have a number of columns equal to 3 (powerpssm, power, and new/old sampling) */
		String str_exp[] = (getValue("sampling_exponents")).split(";"); 
		if (str_exp.length == 0) 
		{
			sampling_exp = new double[1][3];
			sampling_exp[0][0] = 1; /* powerpssm=1 */
			sampling_exp[0][1] = 1; /* power=1 */
			sampling_exp[0][1] = 1; /* 1=new sampling (0=old sampling) */
		}
		
		sampling_exp = new double[str_exp.length][3];
		for(int i = 0; i < str_exp.length; i++) 
		{
			String exps[] = (str_exp[i]).split(" ");
			if (exps.length != 3) 
		    	return "Error: each item of the sampling_exponents variable specified in the \n" +
		    			"configuration file \"" + file_name + "\" must contain exactly 3 values!";				
			
			for(int j = 0; j < 3; j++) 
				try 
				{ 
					sampling_exp[i][j] = Double.parseDouble(exps[j]); 
				}
				catch (Exception e) 
				{
			    	return "Error: one of the values in the sampling_exponent variable (" +  exps[j] +	"),\nas specified in the " +
			    		" configuration file \"" + file_name + "\", is not valid!\nThe phi_prior variable must contain numeric values!";			
				}
	    }
		
		/* *********************************************************************************** */ 
		str = getValue("d");
		try 
		{ 
			d = Double.parseDouble(str); 
		}
		catch (Exception e) 
		{ 
			d = 0; 
		}
		
		if ((d < 0) || (d > 1)) 
			d = 0;
	
		
		str = getValue("pseudocounts_priortype");
		try 
		{ 
			pseudocounts_class = Double.parseDouble(str); 
		}
		catch (Exception e) 
		{ 
			pseudocounts_class = 1; 
		}
		
		if (pseudocounts_class <= 0) 
			pseudocounts_class = 1;
	
		
		str = getValue("pseudocounts_putative_priortype");
		
		try 
		{ 
			pseudocounts_putative_class = Double.parseDouble(str); 
		}
		catch (Exception e) 
		{ 
			pseudocounts_putative_class = 3; 
		}
		
		if (pseudocounts_putative_class <= 0) 
			pseudocounts_putative_class = 3;
		
		str = getValue("flat_prior_other_priortype");
		try 
		{ 
			flat_prior_other_class = Double.parseDouble(str); 
		}
		catch (Exception e) 
		{ 
			flat_prior_other_class = 0.5; 
		}
		
		if ((flat_prior_other_class <= 0) || (flat_prior_other_class >= 1)) 
			flat_prior_other_class = 0.5;
	
		return "";
	}
	
	

	/** Reads the background model from the given file, validates the content
	 * and checks whether or not the order is consistent with the model */
	public static String readBackground()
	{
		/* first store all the lines in an ArrayList */
		BufferedReader br = null;
		ArrayList<String> array = new ArrayList<String>();
		String line;
		
		try 
		{
			br = new BufferedReader(new FileReader(Parameters.back_file));
			
			while ((line = br.readLine()) != null)
			{
				line.trim(); /* add only non-empty lines */
				if (line.length() > 0)
					array.add(line);
			}
			br.close();
		}
		catch (IOException e) 
		{
			try 
			{ 
				br.close();
			} 
			catch (Exception ee) {} /* should not enter here */
			
			return "File error: " + e.getMessage();
		}
		
		/* next we check the size and get the order */
		int max = getMaxOrder(array.size());
		
		if (max < 0)
			return "Error: the number of elements in the background model file \n\"" +
		       Parameters.back_file + "\" (" + array.size() + ") is not valid!\n" +
		       "It must be 4 or 4+4^2 or 4+4^2+4^3 etc.";
	
		
		if (Parameters.bkgrOrder > max) 
			return "Error: the background model order (" + Parameters.bkgrOrder + 
				") is too large! The background file\n\"" + Parameters.back_file + 
				"\" allows a maximum order of " + max + ".";
		
		/* next check: all items must be numbers between 0 and 1 */
		int i = 0;
		back = new double[array.size()];
		
		try 
		{
			for (i=0; i<array.size(); i++) 
			{
				back[i] = Double.parseDouble((String)array.get(i));
				
				if ((back[i] < 0) || (back[i] > 1))
					return "Error: the background model file \"" + Parameters.back_file + 
					"\"\ncontains an invalid value: " + back[i] + 
					"\nAll values should be between 0 and 1!";
			}
		}
		catch (Exception e) 
		{
			return "Error: the background model file \"" + Parameters.back_file + 
			"\"\ncontains an invalid value: \"" + array.get(i) + 
			"\".\nAll values must be real numbers between 0 and 1!"; 			
		}
		
		/* the last check: each group of 4 consecutive numbers should add to 1 */
		double sum;
		
		for (i=0; i<back.length; i+=4) 
		{
			sum = back[i] + back[i+1] + back[i+2] + back[i+3];
			
			if (Math.abs(sum - 1) > precision)
				return "Error: invalid background file \"" + Parameters.back_file + "\"!\n" + 
					"The items from " + i + " to " + (i+3) + " do not sum to 1, but " + sum + "!";
		}
		
		return "";
	}
	
	
	/** Returns the maximum order of a Markov Model, given the size of the table,
	 * or -1 if the size is not valid. */
	public static int getMaxOrder(int bkgrTableSize) 
	{
		if (bkgrTableSize < 4)
			return -1;
		
		int order = 0, sum = 4, prod = 4;
		
		while (sum < bkgrTableSize) 
		{
			order += 1;
			prod *= 4;
			sum += prod;
		}
		
		if (sum == bkgrTableSize)
			return order;
		
		return -1;		
	}
	
	
	/** Reads the names of the TFs from the TF-file and stores them into tf_names. */
	public static String read_TFnames()
	{	
		BufferedReader br = null;
		ArrayList<String> array = new ArrayList<String>(); 
		String line;
		try 
		{
			br = new BufferedReader(new FileReader(Parameters.tf_path));
			
			while ((line = br.readLine()) != null) 
			{
				line.trim(); /* add only non-empty lines */
				if (line.length() > 0)
					array.add(line);
			}
			br.close();
		}
		catch (IOException e) 
		{
			try 
			{ 
				br.close();
			} 
			catch (Exception ee) {} /* should not enter here */
			
			return "File error: " + e.getMessage();
		}
		
		tf_names = null;	
		tf_names = new String[array.size()];
		
		for (int i=0; i<array.size(); i++) 
			tf_names[i] = (String)array.get(i);
	
		return "";		
	}

	
	/**
	 * Sets the output parameters so they become available for vizualisation.
	 * @param Z
	 * @param C
	 * @param phi
	 * @param seq
	 * @param seq_names
	 * @param tf_name
	 */
	public static void setOutput(int Z[], int C[], double phi[][], String seq[], String seq_names[], String tf_name) 
	{
		/* set bestZ */
		bestZ = new int[Z.length];
		for(int i = 0; i < Z.length; i++)
			bestZ[i] = Z[i];

		/* set bestC */
		bestC = new int[C.length];
		for(int i = 0; i < C.length; i++)
			bestC[i] = C[i];
		
		/* set bestPhi */
		bestPhi = new double[phi.length][phi[0].length];
		for(int i = 0; i < phi.length; i++)
			for(int j = 0; j < phi[0].length; j++)
				bestPhi[i][j] = phi[i][j];
	
		/* set comboseq */
		comboseq = new String[seq.length];
		for(int i = 0; i < seq.length; i++)
			if(Parameters.revflag)
				comboseq[i] = GibbsStatic.get_string_from_numbers(seq[i].substring(0, seq[i].length()/2));
			else
				comboseq[i] = seq[i];
		
		/* set comboseq_names */
		comboseq_names = new String[seq_names.length];
		for(int i = 0; i < seq_names.length; i++)
			comboseq_names[i] = seq_names[i];
				
		outputTF_name = tf_name;
		outputParams_are_set = true;
	}
	
	
	/**
	 * 
	 * @return
	 */
	public static String getParameterValuesForCommandLine() 
	{
		String err = Parameters.readBackground();
		if (err.compareTo("") != 0)
			return err;
		
		System.out.println("Background model... OK");
		
		
		/* setting the TF names */
		if (Parameters.individualTF) 
		{ 
			/*we apply the alg for a TF only */
			Parameters.tf_names = new String[1];
			Parameters.tf_names[0] = Parameters.tf_name;
		}
		else 
		{ 
			/* We apply the alg to all the TFs in a file */
			/* the TF names file must be a readable file */
			err = Parameters.read_TFnames();
			if(err.compareTo("") != 0)
				return err;
		}
		
		/* next we have to check that for every TF there is a file in the data path
		 * with the same name as the TF */
		for(int i = 0; i < Parameters.tf_names.length; i++) 
		{
			String path = Parameters.fname_path + "/" + Parameters.tf_names[i] + ".fasta";
			File file = new File(path);
		    
			if (!file.exists()) 
				return "Error: the data file \"" + path + "\"\ncorresponding to the TF \"" + Parameters.tf_names[i] + "\" does not exist!";
			
			if (!file.isFile()) 
				return "Error: the data file \"" + path + "\"\ncorresponding to the TF \"" + Parameters.tf_names[i] + "\" is not a valid file!";
			
			if (!file.canRead()) 
				return "Error: the data file \"" + path + "\"\ncorresponding to the TF \"" + Parameters.tf_names[i] + "\" is not a readable file!";
		}
		
		/* *********************************************************************************** */
		/* SET THE PRIOR INFORMATION */
		
		/* first check single/multiple files */
		if (Parameters.multiple_priors == false) 
			/* uniform or single prior */
			if (Parameters.otherclass == true) 
				Parameters.prior_dirs = new String[0];
		else 
		{ 
			/* multiple priors */
	    	/* the class names are in Parameters.prior_dirs */
	    }
		
	    /* Each entry in Parameters.prior_dirs (if such entries exist) must be the name 
	     * of a readable directory. */	
		for(int i = 0; i < Parameters.prior_dirs.length; i++) 
		{
			File dir = new File(Parameters.prior_dirs[i]);
			if ((!dir.exists()) || (!dir.isDirectory()) || (!dir.canRead()))
				return "Error: the prior directory: \"" + Parameters.prior_dirs[i] + "\" does not exist or is not readable!";
		}

		if (err.compareTo("") != 0)
			return err;
		
		System.out.println("Prior files (" + Parameters.prior_dirs.length + ")... OK");
		return "";
	}
}