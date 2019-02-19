package priority.gibbs;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;
import priority.Parameters;

/**
 * GibbsStatic - the class contains static functions
 * used by the gibbs sampler.
 * @author raluca
 * 
 * Some updates of code have been made by Carlos A. Sierra (carlos.andres.sierra.v@gmail.com)
 */
public class GibbsStatic 
{	
	public static final char DNAchars[] = {'a','c','g','t'};
	public static final double noprior[] = {0,0,0,0};
	
	
	/**
	 * Reads the DNA sequences for a certain TF and returns an array of strings of {0,1,2,3}
	 * @param tf_name
	 * @return
	 * @throws GibbsException
	 */
	public static String[] get_DNAsequences(String tf_name) throws GibbsException
	{
		String fname_file, line, new_line, name_line = null;
		BufferedReader br = null;
		ArrayList<String> local_array = new ArrayList<String>(); 
		
		fname_file = Parameters.fname_path + "/" + tf_name + ".fasta";
		System.out.println("Reading DNA sequences from " + fname_file);
			
		try 
		{
			br = new BufferedReader(new FileReader(fname_file));
			line = br.readLine();

			while (line != null) 
			{
				if (line.length() == 0) 
					line = br.readLine(); 
				else
					if (line.charAt(0) == '>') 
					{
						name_line = line;
						line = br.readLine();
					}
					else 
					{
						/* the next lines should be sequences of acgt -> store as 0123 */
						new_line = "";
						while((line != null) && (line.charAt(0) != '>'))
						{
							for(int i = 0; i < line.length(); i++)
								switch (line.charAt(i)) 
								{
									case 'a': case 'A': new_line = new_line + '0'; break;
									case 'c': case 'C': new_line = new_line + '1'; break;
									case 'g': case 'G': new_line = new_line + '2'; break;
									case 't': case 'T': new_line = new_line + '3'; break;
									default: new_line = new_line + '*';
								}
							
					        line = br.readLine();
						}
						
						local_array.add(name_line);
						local_array.add(new_line);
					}
			}
			
			br.close();
		}
		catch (IOException e) 
		{
			try 
			{ 
				br.close();
			} 
			catch (Exception ee) {}
			
			throw new GibbsException("GibbsException: " + e.getMessage());
		}	
		
		if (local_array.size() < 1)
			throw new GibbsException("GibbsException: no DNA sequences in file \"" + fname_file + "\"!");
		
		String[] string_array = new String[local_array.size()];
		for (int i=0; i<local_array.size(); i++)
			string_array[i] = (String)local_array.get(i);
		
		return string_array;
	}

	
	/**
	 * Reads a positional prior from a fasta-like prior file
	 * @param fname_file
	 * @return
	 * @throws GibbsException
	 */
	public static double[][] get_positional_prior(String fname_file) throws GibbsException
	{
		String line, line2;
		BufferedReader br = null;
		ArrayList<double[]> local_array = new ArrayList<double[]>(); 
		Parameters.wsizeMinFromPrior = -1;
		
		System.out.println("Reading FASTA-like prior from " + fname_file);
			
		try 
		{		
			br = new BufferedReader(new FileReader(fname_file));
			line = br.readLine();
		
			while((line != null) && (line.length() == 0)) 
				line = br.readLine(); 
			
			if (line.startsWith("minmotiflength")) 
			{
				line2 = br.readLine();
			
				if (line2.startsWith("maxmotiflength")) 
				{
					String[] result = line.split(" |=");
					String[] result2 = line2.split(" |=");
					int min = -1; int max = -1;
				
					try
					{
						min = Integer.parseInt(result[ result.length - 1 ]);
						max = Integer.parseInt(result2[ result2.length - 1 ]);
					} 
					catch (Exception e) 
					{ 
						System.out.println("Error in a prior file: " + e + " (line ignored)"); 
				    }
					
					if (min >= Parameters.wsizeMin && min <= Parameters.wsizeMax &&	max >= Parameters.wsizeMin && max <= Parameters.wsizeMax && min <= max) 
					{
						/* the min and max values specified in the prior file are valid */
						Parameters.wsizeMinFromPrior = min;
						Parameters.wsizeMaxFromPrior = max;
					}
					
					line = br.readLine();
				}
			}
					
			while(line != null) 
			{
				if(line.length() != 0)
					if (line.charAt(0) == '>') {}
					else
					{
						/* this line should be a sequence of probabilities -> store it */
						String[] result = line.split(" ");
						double array[] = new double[result.length];
						
						try
						{
							for(int i = 0; i < result.length; i++)
								array[i] = Double.parseDouble(result[i]);
						} 
						catch (Exception e) 
						{
							e.printStackTrace(); 
						}
						
						local_array.add(array);
					}
				
				line = br.readLine();
			}			
			br.close();
		}
		catch (IOException e) 
		{
			try { 
				br.close();
			} 
			catch (Exception ee) {}
			throw new GibbsException("GibbsException: " + e.getMessage());
		}	
		
		if (local_array.size() < 1)
		{
			throw new GibbsException("GibbsException: no FASTA-like prior in file \"" +	fname_file + "\"!");
		}
			
		double final_array[][] = new double[local_array.size()][];
		
		for(int i = 0; i < local_array.size(); i++)
			final_array[i] = (double[])local_array.get(i);
		
		return final_array;
	}
	
	
	/**
	 * Transforms a string of {0,1,2,3} into a string og {a,c,g,t}.
	 * @param seq
	 * @return
	 */
	static public String get_string_from_numbers(String seq) 
	{
		String new_seq = "";
		
		for(int i = 0; i < seq.length(); i++)
			switch (seq.charAt(i)) 
			{
				case '0': new_seq = new_seq + 'A'; break;
				case '1': new_seq = new_seq + 'C'; break;
				case '2': new_seq = new_seq + 'G'; break;
				case '3': new_seq = new_seq + 'T'; break;
				case '*': new_seq = new_seq + 'N'; break;
			}
						
		return new_seq;
	}
	
	
	/**
	 * Returns the reverse complement of a seq (acgt = 0123)
	 * @param seq
	 * @return
	 */
	static public String get_reverse(String seq)
	{
		String new_seq = "";
		
		for(int i = 0; i < seq.length(); i++)
			switch(seq.charAt(i))
			{
				case '0': new_seq = '3' + new_seq; break;
				case '1': new_seq = '2' + new_seq; break;
				case '2': new_seq = '1' + new_seq; break;
				case '3': new_seq = '0' + new_seq; break;
				case '*': new_seq = '*' + new_seq; break;
			}
		
		return new_seq;
	}
	
	
	/**
	 * Computes the index, in the precomputed table representing a class prior, of a string of {0,1,2,3}
	 * @param seq
	 * @return
	 */
	static public int get_index(String seq)
	{
		int index = 0;
		
		for (int i=0; i<seq.length(); i++)
			index = index * 4 + (int)(seq.charAt(i) - '0');
				
		return index;
	}
	
	
	/**
	 * Calculates phi (from the pseudo counts + the counts of the sites contributing to the current alignment Z[-index]
	 * @param Z
	 * @param seq
	 * @param index
	 * @param wsize
	 * @param qprior
	 * @return
	 */
	static public double[][] calPhi(int Z[], String seq[], int index, int wsize, double qprior[])
	{
		double[][] phi = new double[4][wsize];
		double total = 0;
		int i, j;
		
		/* initialize phi with the pseudocounts */
		for(i = 0; i < 4; i++) 
		{
			total += qprior[i];
			
			for(j = 0; j < wsize; j++)
				phi[i][j] = qprior[i];
		}		
		
		/* count the occurrences for each nucleotide on each position */ 
		for(i = 0; i < seq.length; i++)
			if(i != index)  /* except the seq used in the current iteration */
				if (Z[i] >= 0) 
				{ /* if there is an occurrence of a motif in this sequence */
					for(j = 0; j < wsize; j++)
						phi[seq[i].charAt(Z[i]+j) - '0'][j] += 1;
					
					total++;
				}
		
		
		/* normalize */
		for(i = 0; i < 4; i++) 
			for(j = 0; j < wsize; j++)
				phi[i][j] /= total;		
				
		return phi;
	}
	
	
	/**
	 * Calculates:  P(array|phi) / P(array|background)
	 * @param phi
	 * @param wsize
	 * @param array
	 * @param back
	 * @param order
	 * @return
	 */
	static public double calA(double phi[][], int wsize, String array, double back[], int order)
	{
		int index, j, cnt, offset;
		double A = 1;
	
		/* the prob. of the subseq according to the PSSM */
		for(j = 0; j < wsize; j++)
			A *= phi[array.charAt(j)-'0'][j];
		
		/* the prob. of the subseq | background */ 
		for(j = 0; j < wsize; j++)
		{
			index = 0;
			cnt = 0;
			offset = 0;
			
			for(int k = Math.max(0,j-order); k <= j; k++)
			{
				index = index * 4 + (array.charAt(k) - '0');
				offset += 1 << (2 * cnt); // 4^cnt = 2^(2*cnt);
				cnt++;
			}
			
			index = index + offset - 1;
			A /= back[index];
		}	
		
		return A;
	}
	
	
	/**
	 * Calculates:  P(array|phi)
	 * @param phi
	 * @param wsize
	 * @param array
	 * @param back
	 * @param order
	 * @return
	 */
	static public double calonlypssm(double phi[][], int wsize, String array, double back[], int order)
	{
		double A = 1;
	
		// the prob. of the subseq according to the PSSM 
		for(int j = 0; j < wsize; j++)
			A *= phi[array.charAt(j)-'0'][j];
				
		return A;
	}
	
	
	/**
	 * This function is used for computing the log likelihood.
	 * @param phi
	 * @param phi_prior
	 * @param gamma
	 * @param cprior
	 * @return
	 */
	public static double logl_phi_gamma(double phi[][], double phi_prior[], double gamma[], double cprior[])
	{
		double prod, logl = 0;
		
		// the phi part 
		for(int i = 0; i < 4; i++)
		{
			prod = 1;
			
			for (int j=0; j<phi[0].length; j++)
				prod *= phi[i][j];
			
			logl += (phi_prior[i] - 1) * Math.log(prod);
		}
		
		// the gamma part 
		for(int cl = 0; cl < gamma.length; cl++)
			logl += (cprior[cl] - 1) * Math.log( gamma[cl] );
		
		return logl;
	}


	/**
	 *   
	 * @param cr_class = the current class
	 * @param nc = the number of classes
	 * @param cprior = the pseudo counts for each class
	 * @param C = the class assigned to each sequence
	 * @param index = the index of the current sequence
	 * @return
	 */
	public static double prob_gen_prior(int cr_class, int nc, double cprior[], int C[], int index)
	{
		double[] prob = new double[nc];
		double sum = 0;
		int i = 0;
		
		// initialize with the pseudo counts 
		for(i = 0; i < nc; i++)
			prob[i] = cprior[i];	
		
		// add the actual counts 
		for(i = 0; i < C.length; i++)
			if (index != i && C[i] > -1)
				prob[ C[i] ] += 1;
		
		for(i = 0; i < nc; i++)
			sum += prob[i];
		
		return prob[cr_class] / sum;
	}

	
	/**
	 * Randomly generates an index from 0 to W.length-1, according to the weights W and their sum sumW
	 * @param W
	 * @param sumW
	 * @return
	 */
	public static int rand_sample(double W[], double sumW)
	{
		/* generate a number between 0 and sumW */
		if (sumW == 0)
			return -1;
		
		double value = (new Random()).nextDouble() * sumW, sum = 0.0;
		int index = 0;
		
		while (sum <= value && index < W.length)
			sum += W[index++];
		
		return index-1;
	}
	
	
	/**
	 * Returns the max value in the array
	 * @param array
	 * @return max_value
	 */
	public static double max(double array[])
	{
		double x = array[0];
		
		for(int i = 1; i < array.length; i++)
			if (array[i] > x)
				x = array[i];
		
		return x;
	}

	
	
	/* ******************** Functions used for printing ******************** */
	
	/**
	 * Computes the number of occurrences for each class.
	 * @param C
	 * @param Z
	 * @param nc
	 * @return
	 */
	public static int[] class_counts(int[] C, int[] Z, int nc)
	{
		int[] counts = new int[nc+1];
		int i;
		
		for(i = 0; i < nc + 1; i++)
			counts[i] = 0;
		
		for(i = 0; i < C.length; i++) 
			if(Z[i] >= 0)
				counts[C[i]]++;  /* occurrence of a motif for a TF of class C[i] */ 
		
		for(i = 0; i < Z.length; i++)
			if(Z[i] < 0) /* no occurrence */
				counts[nc]++;
		
		return counts;
	}
	
	
	/**
	 * Return a string of exactly <times> characters <c>.
	 * @param c
	 * @param times
	 * @return char repetitions
	 */
	public static String repeatChar(char c, int times)
	{
		String str = "";
		
		for(int i = 0; i < times; i++)
			str = str + c;
		
		return str;
	}
	
	
	/**
	 * Return a formatted version of a double between 0 and 1.
	 * @param x
	 * @param decimals
	 * @return formatted number
	 */
	public static String formatDouble01(double x, int decimals)  
	{
		String str = "";
		long tmp = Math.round(x * Math.pow(10, decimals));
	
		for(int i = 0; i < decimals; i++) 
		{
			str = (tmp % 10) + str;
			tmp = tmp / 10;
		}
		
		return tmp + "." + str;
	}
	
	
	/**
	 * Returns a formatted version of an integer of max 2 digits.
	 * @param x
	 * @param positions
	 * @return formatted number
	 */
	public static String formatInt(int x, int positions) 
	{
		int pos = x < 10 ? positions - 1 : positions - 2;
		
		return repeatChar(' ', pos - pos / 2) + x + repeatChar(' ', pos / 2); 
	}
}