package priority;
import priority.gui.MainWindow;
import priority.gibbs.GibbsRun;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.Color;
import javax.swing.JOptionPane;
import javax.swing.JFrame;
import javax.swing.JWindow;
import javax.swing.ImageIcon;
import javax.swing.JLabel;

/**
 * Priority - the main class
 * It runs both the GUI and the algorithm. 
 * @author raluca
 */
public class Priority 
{
	public MainWindow win;
	public GibbsRun gibbs;
    static public boolean useInterface = true;
    static public String userConfigFile = null;
	static public Image icon = null;
	
	/** Number of milliseconds the splash screen will stay visible. */
	private static int SPLASH_VIEW_TIME = 800;
	
	/** Number of milliseconds that will elapse before beginning the
	 * actual loading of the graphical interface  */
	private static final int WAIT_BEFORE_START = 600;
		
	/** The main object */
	static Priority priority;
	
	/** The main function. */
	public static void main (String args[]) throws Exception 
	{	
		/* Create the main object. */
		priority = new Priority();
		
		if (args.length == 0) { /* continue (run the GUI version) */		
		} 
		else if (args.length == 1) {
			if (args[0].equals("-nogui"))
				useInterface = false;
			else
				printHelp();
		}
		else if (args.length == 2) {
			/* the parameters should be "-f" and the name of a config file */
			if (args[0].equals("-f")) {
				userConfigFile = args[1];
			}
			else
				printHelp();
		}
		else if (args.length == 3) {
			/* the parameters should be "-nogui", "-f" and the name of a config file */
			if (!args[0].equals("-nogui"))
				printHelp();
			else {
				useInterface = false;	
				if (args[1].equals("-f")) {
					userConfigFile = args[2];
				}
				else
					printHelp();
			}
		} 
		else
			printHelp();
		
		if (useInterface) {
			/* Set the icon. */
			String imgName = "/priorityicon.jpg";
			java.net.URL url = priority.getClass().getResource(imgName);
			if (url != null) {
				ImageIcon picture = new ImageIcon(url);
				icon = picture.getImage();
			}

			/* This thread loads the splash screen, allowing the application to load in the background */
			new Thread() {
				public void run() {
					JWindow splash = new JWindow();                    
	
					String imgName = "/prioritysplash.jpg";
					java.net.URL url = this.getClass().getResource(imgName);
					ImageIcon picture = new ImageIcon(url);
	
					JLabel label = new JLabel(picture);
					int w = picture.getIconWidth();
					int h = picture.getIconHeight();
					splash.setSize(w, h);
					splash.getContentPane().add(label, BorderLayout.CENTER);
					Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
					splash.setLocation(
							(dim.width - w) / 2,
							(dim.height - h) / 2);
					splash.setVisible(true);
					
					try {
						Thread.sleep(SPLASH_VIEW_TIME);
					} catch (InterruptedException e) {}					
					splash.dispose();				
				}
			}.start();	
		}
		
		new Thread() {
			public void run() {
				if (useInterface) {
					try {
						Thread.sleep(WAIT_BEFORE_START);
					} catch (InterruptedException e) {}
				}
				
				/* Read the strings for the GUI. */
				String err = Strings.setStrings();
				if (err.compareTo("") != 0) {
					System.out.println(err);
					if (useInterface) {
						JOptionPane.showMessageDialog(null, err, "Configuration error", JOptionPane.ERROR_MESSAGE);
					}
					System.exit(0);
				}

				/* Read and set the default parameters. */
				if (userConfigFile == null)
					err = Parameters.setDefaultParameters();
				else
					err = Parameters.setDefaultParameters(Priority.userConfigFile);
				if (err.compareTo("") != 0) {
					System.out.println(err);
					if (useInterface) {
						JOptionPane.showMessageDialog(null, err, "Config error", JOptionPane.ERROR_MESSAGE);
					}
					System.exit(0);
				}
				
		        if (useInterface) {
					/*Makes the JFrame.setMinimumSize work */
					try {
						String os = "unknown";
						try { os = System.getProperty("os.name");
						} catch (SecurityException e) {}
						if (os.toLowerCase().indexOf("linux") < 0)
							javax.swing.UIManager.setLookAndFeel(javax.swing.UIManager.getSystemLookAndFeelClassName());
						else
							JFrame.setDefaultLookAndFeelDecorated(true);
			        } catch(Exception e) {}

					/* Create and show the main window. */
					priority.win = new MainWindow(priority);
					priority.win.init();
					priority.win.setVisible(true);
				}
				else {
					Parameters.getParameterValuesForCommandLine();
					priority.startGibbs();
				}
			}
		}.start();
	}
	
	/** Starts a thread that runs the gibbs sampler. */
	public void startGibbs() {
		gibbs = new GibbsRun(this);
		gibbs.set_local_params();
		gibbs.start();	
	}

	/** Stops the thread created in startGibbs. */
	public void stopGibbs() {
		if ((gibbs != null) && gibbs.isAlive())
			gibbs.stopThread();
		gibbs = null;
		System.gc();
	}

	/** This function is called by the gibbs sampler when it is done. */  
	public void finished() {
		win.activateStart();
	}
	
	/** Prints the output. */
	public void printOutput(String str) {
		this.win.appendTextToTextPanel(str);
	}
	
	/** Prints the output using a specific color. */
	public void printOutput(String str, Color color) {
		this.win.appendTextToTextPanel(str, color);
	}
	
	public static void printHelp() {
		System.out.println("Usage:  ");
		System.out.println("  java -jar priority.jar                         Run PRIORITY (GUI version)");
		System.out.println("  java -jar priority.jar -nogui                  Run PRIORITY (command line version)");
		System.out.println("  java -jar priority.jar -f config_file          Run PRIORITY (GUI) with the specified config file");
		System.out.println("  java -jar priority.jar -nogui -f config_file   Run PRIORITY (command line) with the specified config file");
		System.out.println("  java -jar priority.jar -h                      Print this help message.");
		System.exit(0);
	}
}
