package priority.gibbs;

/**
 * GibbsException - the exceptions thrown from the 
 * gibbs sampler will be of type GibbsException.
 * @author raluca
 */
public class GibbsException extends Exception 
{
	private static final long serialVersionUID = 1;
	public GibbsException() {}
	
	public GibbsException(String msg) 
	{
		super(msg);
	}
}