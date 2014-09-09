package priority.gibbs;

import priority.gibbs.GibbsStatic;
import java.util.Random;

/**
 * This class represents the behavior of a Trial, based on threads.
 * @author Eng. Paula Reyes, Ph.D. - M.Sc. Eng. Carlos Sierra 
 * Univerisdad Antonio Narino - Colombia
 */
public class Trial extends Thread {

	/* Running parameters */
	private double phi_prior[];      /* dirichlet prior counts for the model parameters (phi) */
	private double[][] sampling_exp; /* the exponents for prior and likeliood, used for sampling */
	private double[] back;           /* the background model */
	
	private int wsize;                  /* the window size */	
	private int nc;                     /* the number of classes (prior_dirs.length + 1 for the "other" class) */	
	private double comboprior[][][];    /* the priors for each class/prior, each sequence, each position (different for each TF) */
	private String comboseq[];          /* the sequences for the current TF (as strings of {0,1,2,3} <-> {A,C,G,T} */
	
	/* Other variables */
	private double cprior[];     /* dirichlet prior counts for the classes */
	private double denom[][];    /* the normalization constant (for each class and each sequence) */
	
	/* The variables we sample */
	public int Z[], bestZ[] = null;
	public int C[], bestC[] = null;
	
	/* More variables */
	private boolean noocflag;
	private int outputStep; 
	private int bkgrOrder;
	private int iter; 
	public int tr;
	private int trials;
	public double[] logl;
	public double bestlogl;
	public int bestlogl_iter;
	
	
	
	/**
	 * Constructor of class. Receives all data to run an independent trial.
	 */
	public Trial(String[] comboseq, double[][][] comboprior, int nc, int trials, double[] cprior, double[] phi_prior, double[] back, 
			int iter, double[][] sampling_exp, boolean noocflag, int outputStep, int bkgrOrder, int tr, int[] Z, double[][] denom, int wsize, int tf)
	{
		this.comboseq = comboseq; 
		this.comboprior = comboprior; 
		this.nc = nc; 
		this.trials = trials; 
		this.cprior = cprior; 
		this.phi_prior = phi_prior; 
		this.back = back; 
		this.iter = iter; 
		this.sampling_exp = sampling_exp; 
		this.noocflag = noocflag; 
		this.outputStep = outputStep; 
		this.bkgrOrder = bkgrOrder; 
		this.tr = tr; 
		this.Z = Z; 
		this.denom = denom; 
		this.wsize = wsize; 
		
		bestlogl = -Double.MAX_VALUE;
		bestlogl_iter = -1;				
		logl = new double[iter];
	}
	
	
	/**
	 *  This method runs all steps for Trial, including iterations, sampling, and selection of best solution found.
	 */
	private void execute()
	{
		double[][] phi;
		Z = new int[comboseq.length];
		bestZ = new int[comboseq.length];
		C = new int[comboseq.length];
		bestC = new int[comboseq.length];
			
		int[] C = new int[comboseq.length]; 
			
		Random rand = new Random();
		
		/* initialize C and Z (uniform probability) */
		for (int i=0; i<comboseq.length; i++)
		{
			C[i] = rand.nextInt(nc); /* random integers between 0 and nc-1 */
			
			/* generate a weight vector that has 0 wherever the prior for the
			 * class C[i] has 0, and 1 everywhere else */
			double temp[] = new double[comboseq[i].length()];
			double tempsum = 0;
			
			for (int j=0; j<comboseq[i].length(); j++) 
			{
				temp[j] = comboprior[ C[i] ][i][j];
				
				if (temp[j] > 0)
				{
					temp[j] = 1;
				}
				
				tempsum += temp[j];
			}
			
			Z[i] = GibbsStatic.rand_sample(temp, tempsum);
			if (Z[i] == -1) /* debugging only */
			{
				System.out.println("ATTN: rand_sample returned -1!!!");
			}
		}
				
		int index = rand.nextInt(comboseq.length); /* choose random index */
		
		/* run with different values for power and powerpssm */
		int power = 1;
		int powerpssm = 1;
		boolean flag_for_oldsampling = false;
		
		int n_sampling = this.sampling_exp.length;
		
		for (int i_sampling=1; i_sampling <= n_sampling; i_sampling++) 
		{
			if (tr < i_sampling * trials / n_sampling) 
			{
				powerpssm = (int)sampling_exp[i_sampling-1][0];
				power = (int)sampling_exp[i_sampling-1][1];
			
				if ((int)sampling_exp[i_sampling-1][2] == 0)
				{
					flag_for_oldsampling = true;
				}
				
				break;
			}
		}				
		
		if (powerpssm != 1 || power != 1) 
		{
			noocflag = false;
		}
			
		for (int cnt=0; cnt<iter; cnt++)
		{
			if((cnt+1)%outputStep == 0)
			{
				power=Math.max(1,power-3);
				powerpssm=Math.max(1,powerpssm-1);
		
				if(power > 1 || powerpssm >1)
				{
					noocflag = false;
				}
				else
				{
					noocflag = true;
				}
			}
			
			/* pick the next sequence (+1 or randomly: index = rand.nextInt(comboseq.length))*/
			index = (index+1) % (comboseq.length);
			
			/* calculate the current estimate of phi */
			phi = GibbsStatic.calPhi(Z, comboseq, index, wsize, phi_prior);
			
			/* next we have to sample Z[index], but before that we must compute the weights:*/
			int length = comboseq[index].length();
			double[] W = new double[length+1];
			double sumW = 0;
			
			/* we initialize the weight vector for the class also (will be used later) */
			double[] Wc = new double[nc];
			double sumWc = 0;
			
			for (int j=0; (j<length); j++)  // && (stop_thread==false); j++) 
			{	
				if (comboprior[C[index]][index][j] == 0) 
				{
					W[j] = 0;
				}
				else 
				{ 
					double quantity1 = GibbsStatic.calA(phi, wsize, 
							comboseq[index].substring(j,j+wsize), back, bkgrOrder);
					
					if(flag_for_oldsampling) 
					{
						quantity1=quantity1*Math.pow(GibbsStatic.calonlypssm(phi, wsize, 
							comboseq[index].substring(j,j+wsize), back, bkgrOrder),powerpssm-1);
					}
					else 
					{
						quantity1=Math.pow(quantity1,powerpssm);
					}
					
					double quantity2 = Math.pow(comboprior[C[index]][index][j],power) / (1-Math.pow(comboprior[C[index]][index][j],power));
					W[j] =  quantity1 * quantity2;							
					sumW += W[j];
				}
			}
								
			if (GibbsStatic.max(W) <= 0)
			{
				Z[index] = -1; /* debugging case!!! "all weights 0") Should not happen!!! */
			}
			else 
			{
				if (noocflag)
				{
					W[length] = 1; /* to sample the "no motif" case */
				}
				else
				{
					W[length] = 0;
				}
				
				sumW += W[length];
				
				/* now we finally sample Z[index] */
				Z[index] = GibbsStatic.rand_sample(W, sumW);
				
				if(Z[index]==-1)
				{
					System.out.println("ATTN "+index+" "+Z[index]+" "+length);
				}
				
				/* next we compute the weights for sampling C[index] */										
				if (Z[index] == length) /* the sampler picked the appended 1 ("no motif") */ 
				{
					Z[index] = -1; 
					
					for (int c=0; c<nc; c++) 
					{
						Wc[c] = (1/denom[c][index]) * GibbsStatic.prob_gen_prior(c,nc,cprior,C,index);
					}
				} 
				else /* we set the weights to sample class where a motif exists */ 
				{
					for (int c=0; c<nc; c++)
					{
						if(Z[index]>length-wsize+1 || Z[index]<0)
						{
							System.out.println(index+" "+Z[index]+" "+length);
						}
						
						Wc[c] = comboprior[c][index][Z[index]] * GibbsStatic.prob_gen_prior(c,nc,cprior,C,index);
					}
				}
				
				sumWc = 0;
				
				for (int c=0; c<nc; c++)
				{
					sumWc += Wc[c];
				}
				
			} /* end if (GibbsStatic.max(W) == 0) - the debug case */
			
			/* now we sample C[index] */
			C[index] = GibbsStatic.rand_sample(Wc, sumWc);
			
			/* recompute phi */
			phi = GibbsStatic.calPhi(Z, comboseq, -1, wsize, phi_prior);
			
			/* now the MAP calculation: */
			double gamma[] = new double[nc];
			double sum_gamma = 0;
			
			for (int cl=0; cl<nc; cl++) 
			{
				gamma[cl] = GibbsStatic.prob_gen_prior(cl,nc,cprior,C,-1);
				sum_gamma += gamma[cl];
			}
			
			/* normalize it => a prob distribution */
			for (int cl=0; cl<nc; cl++)
			{
				gamma[cl] = gamma[cl] / sum_gamma;
			}
			
			
			/* sum log(P(gamma)) and log(P(phi)) */
			logl[cnt] = GibbsStatic.logl_phi_gamma(phi, phi_prior, gamma, cprior);
			
			for (int i=0; i<Z.length; i++)
			{	
				if (Z[i] > -1) /* there is an occurrence of the motif */
					logl[cnt] = logl[cnt] +
					   Math.log(GibbsStatic.calA(phi, wsize, comboseq[i].substring(Z[i],
							   Z[i]+wsize), back, bkgrOrder)) +
					   Math.log(comboprior[C[i]][i][Z[i]]) - 
					   Math.log(1-comboprior[C[i]][i][Z[i]]) - 
					   Math.log(denom[C[i]][i]) + 
					   Math.log(GibbsStatic.prob_gen_prior(C[i], nc, cprior, C, -1));
				else /* no occurrence */
					logl[cnt] = logl[cnt] -
					   Math.log(denom[C[i]][i]) +
					   Math.log(GibbsStatic.prob_gen_prior(C[i], nc, cprior, C, -1));
			}
			
			
			if (logl[cnt] > bestlogl) 
			{
				bestlogl = logl[cnt];
				bestlogl_iter = cnt;
				
				for (int i=0; i<Z.length; i++) 
				{
					bestZ[i] = Z[i];
				}
				
				for (int i=0; i<C.length; i++) 
				{
					bestC[i] = C[i];
				}
			}
		} /* end for (int cnt=0; cnt<iter; cnt++) */
	}

	
	/**
	 * This method starts the Trial thread
	 */
	public void run()
	{
		//Call to function that execute Trial's workflow
		this.execute();
	}
	
}