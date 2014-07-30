package priority.gui;
import priority.Priority;
import priority.Strings;
import priority.Parameters;
import priority.gui.view.MotifView;
import priority.gui.view.PositionsView;
import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Dimension;
import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JLabel;
import javax.swing.JButton;
import javax.swing.BorderFactory;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JEditorPane;
import javax.swing.WindowConstants;
import javax.swing.ImageIcon;
import java.io.IOException;

/**
 * The main window of the application.
 * @author raluca
 * (some code borrowed from Matt Edwards, Jason Bosko)
 */
public class MainWindow extends JFrame implements ActionListener 
{
	static final long serialVersionUID = 1;
	
	JPanel contentpane;
	JPanel center, south, west;
	JPanel buttons;
	JButton start, stop, cleartext, exit;
	
	ParamPanel paramPanel;
	TextPanel textPanel;

	Priority mainApp;
	
	
	/** Constructor */
	public MainWindow(Priority priority)
	{
		super(Strings.getString("mainWindowTitle"));	
		mainApp = priority;
	}
	
	
	/** Appends the text to the text panel. */
	public void appendTextToTextPanel(String text) 
	{
		textPanel.appendText(text);
	}
	
	
	/** Appends the text to the text panel (using the specified color. */
	public void appendTextToTextPanel(String text, Color color) 
	{
		textPanel.appendText(text, color);
	}
	
	
	/** Creates and initilizes the components. */
	public void init() 
	{
		this.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				mainApp.stopGibbs();
				System.exit(0);
			}
		});	
		
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		//this.setPreferredSize(new Dimension(900, 750));
		this.setMinimumSize(new Dimension(600, 630));
		//this.setMaximumSize(new Dimension(900, 900)); /* doesn't work */
		
		this.setResizable(true);
		
		//this.setSize(new Dimension(900, 640));
		this.setIconImage(Priority.icon);
		 
		contentpane = (JPanel)this.getContentPane();
		contentpane.setLayout( new BorderLayout() );
		contentpane.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
		
		center = new JPanel( new BorderLayout() );
		center.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createEmptyBorder(6,0,6,6),
				BorderFactory.createLoweredBevelBorder())
				);
		contentpane.add(center, BorderLayout.CENTER);	
		
		south = new JPanel( new BorderLayout() );
		south.setBorder(BorderFactory.createEmptyBorder(0,5,5,5));
		south.setPreferredSize( new Dimension(830, 30));
		contentpane.add( south, BorderLayout.SOUTH);

		west = new JPanel( new BorderLayout() );
		west.setBorder(BorderFactory.createEmptyBorder(1,5,5,5));
		contentpane.add( west, BorderLayout.WEST);
		
		paramPanel = new ParamPanel();
		paramPanel.init();	
		west.add(paramPanel,BorderLayout.CENTER);
		
		textPanel = new TextPanel();
		textPanel.init();
		
		center.add(textPanel, BorderLayout.CENTER);
			
        buttons = new JPanel( new GridLayout(1,3, 5, 0));
        start = new JButton(Strings.getString("start"));
        setButton(start, buttons);		
        stop = new JButton(Strings.getString("stop"));
        setButton(stop, buttons);
        cleartext = new JButton(Strings.getString("cleartext"));
        setButton(cleartext, buttons);
        exit = new JButton(Strings.getString("exit"));
        setButton(exit, buttons);
		south.add(buttons, BorderLayout.CENTER);
		
		MenuMaker.makeMenus(this);
		MenuMaker.viewMenuViewPositions.setEnabled(true);
		MenuMaker.viewMenuViewPSSM.setEnabled(true);
		
		if (Parameters.multiple_priors == true)
		{
			MenuMaker.priorsMenuMultiplePriors.setSelected(true);
		}
		else
		{
			MenuMaker.priorsMenuSinglePrior.setSelected(true);
		}
		
		pack();		 
    } 		

	
	/** Sets button properties. */
	private void setButton(JButton button, JPanel panelButtons) 
	{
		String os = "unknown";
		try 
		{ 
			os = System.getProperty("os.name");
		} 
		catch (SecurityException e) {}
		
		if(os.indexOf("Mac") >= 0) 
		{
			button.setBorder(BorderFactory.createRaisedBevelBorder());
		}
		
		button.addActionListener(this);	
		panelButtons.add(button);
	}

	
	/** Handles the button clicks. */
	public void actionPerformed(ActionEvent ev) 
	{
		if (ev.getSource() == start) 
		{
			String err = paramPanel.setParameterValues();
			if (err.compareTo("") == 0) 
			{
				deactivateStart();
				Parameters.outputParams_are_set = false;
				mainApp.startGibbs();
			}
			else 
			{		
				JOptionPane.showMessageDialog(this, err, "Parameter error", JOptionPane.ERROR_MESSAGE);
			}
 			return;
 		}
 		
		if (ev.getSource() == stop) 
		{
 			mainApp.stopGibbs();
 			activateStart();
 			return;
 		}
		
		if (ev.getSource() == cleartext) 
		{
			this.textPanel.resetText();
 			return;
 		}
		
		if (ev.getSource() == exit) 
		{
			closeApplication();
 			return;
 		}
	}
	
	
	/** Activates the start button */
	public void activateStart() 
	{
		start.setEnabled(true);
		MenuMaker.viewMenuViewPositions.setEnabled(true);
		MenuMaker.viewMenuViewPSSM.setEnabled(true);
		MenuMaker.priorsMenuMultiplePriors.setEnabled(true);
		MenuMaker.priorsMenuSinglePrior.setEnabled(true);		
	}
	
	
	/** Deactivates the start button */
	public void deactivateStart() 
	{
		start.setEnabled(false);
		MenuMaker.viewMenuViewPositions.setEnabled(false);
		MenuMaker.viewMenuViewPSSM.setEnabled(false);
		MenuMaker.priorsMenuMultiplePriors.setEnabled(false);
		MenuMaker.priorsMenuSinglePrior.setEnabled(false);
	}

	
	/** Closes the application.*/
	public void closeApplication() {
		mainApp.stopGibbs();
		this.dispose();
		System.exit(0);
	}
	
	
	/** Saves the current parameters into a file using serialization.*/
	public void saveParametersFile(String fileName) 
	{
		System.out.println("Save all params into file: "+fileName);
	}
	
	
	/** Loads all the parameters from a file using serialization.*/
	public void loadParametersFile(String fileName) 
	{
		System.out.println("Load all params from file: " + fileName);
	}
	
		
	/** Switches to Uniform/SinglePrior view. */
	public void switchToSinglePrior() 
	{
		if (Parameters.multiple_priors == true) 
		{
			this.textPanel.resetText();
			this.paramPanel.switchToSinglePrior();
			this.pack();
			this.validate();
		}
	}

	
	/** Switches to MultiplePriors view. */
	public void switchToMultiplePriors() 
	{
		if (Parameters.multiple_priors == false) 
		{
			this.textPanel.resetText();
			this.paramPanel.switchToMultiplePriors();
			this.pack();
			this.validate();
		}
	}
	
	
	/** Shows the final motif. */
	public void viewFinalPSSM() 
	{
		if (Parameters.outputParams_are_set == false) 
		{
			JOptionPane.showMessageDialog(this, 
					"Please run the algorithm for at least one trial to compute a motif!", 
					"Error", JOptionPane.ERROR_MESSAGE);
			return;
		}
		
		MotifView motifFrame = new MotifView("PRIORITY: Best scoring motif (TF: " + Parameters.outputTF_name + ")");
		motifFrame.init(Parameters.bestPhi);
		
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		motifFrame.validate();
		motifFrame.setLocation((dim.width - motifFrame.getWidth()) / 2, (dim.height - motifFrame.getHeight()) / 2);
		motifFrame.setVisible(true);	
	}
	
	
	/** Shows the occurrences of the final motif. */
	public void viewFinalMotifPositions() 
	{
		if (Parameters.outputParams_are_set == false) 
		{
			JOptionPane.showMessageDialog(this, 
					"Please run the algorithm for at least one trial to compute a motif!", 
					"Error", JOptionPane.ERROR_MESSAGE);
			return;
		}
		
		PositionsView posFrame = new PositionsView("PRIORITY: Motif occurrences (TF: " + Parameters.outputTF_name + ")");
		posFrame.init(Parameters.bestZ, Parameters.comboseq, Parameters.comboseq_names);
		
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		posFrame.validate();
		posFrame.setLocation((dim.width - posFrame.getWidth()) / 2, (dim.height - posFrame.getHeight()) / 2);
		posFrame.setVisible(true);		
	}
	
	
	/** Will show help topics in the future. */
	public void viewHelpTopics() 
	{
		System.out.println("show help topics");
	}
	
	
	/** Shows the file README.txt. */
	public void viewHelpReadme() 
	{
		JFrame readme = new JFrame("PRIORITY: README.txt");
		readme.setSize(630, 500);
		readme.setIconImage(Priority.icon);
		
		
		String fileName = "/README.txt";
		java.net.URL url = this.getClass().getResource(fileName);
	
		if (url != null) 
		{
			JEditorPane readmeArea = null;
			try 
			{
				readmeArea = new JEditorPane(url);
			}
			catch (IOException e) 
			{
				JLabel l = new JLabel("Could not load README.txt file.");
				l.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
				readme.getContentPane().add(l);
				readme.setSize(280, 120);		
			}
			
			readmeArea.setMargin(new Insets(5,5,5,5));
			readmeArea.setCaretPosition(0);
			readmeArea.setContentType("text/plain");
			readmeArea.setEditable(false);
			Font font = readmeArea.getFont();
			readmeArea.setFont(new Font("MonoSpaced", Font.PLAIN, font.getSize()));
			readme.getContentPane().add(new JScrollPane(readmeArea));

		} 
		else 
		{
			JLabel l = new JLabel("Could not load README.txt file.");
			l.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
			readme.getContentPane().add(l);
			readme.setSize(280, 120);
		}
		
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		readme.validate();
		readme.setLocation((dim.width - readme.getWidth()) / 2, (dim.height - readme.getHeight()) / 2);
		readme.setVisible(true);
	}
	
	
	/** Shows the "about" information. */
	public void viewHelpAbout() 
	{
		JFrame about = new JFrame("PRIORITY: About");
		about.setResizable(false);
		about.setIconImage(Priority.icon);
		about.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		about.getContentPane().setBackground(Color.WHITE);
		
		String imgName = "/priorityabout.jpg";
		java.net.URL url = this.getClass().getResource(imgName);
		ImageIcon picture = new ImageIcon(url);
		JLabel label = new JLabel(picture);
		
		Dimension dim1 = Toolkit.getDefaultToolkit().getScreenSize();
		
		int w = picture.getIconWidth() + 4;
		int h = picture.getIconHeight() + 25;
		
		about.setSize(w, h);
		about.getContentPane().add(label, BorderLayout.CENTER);
		about.setLocation(
				(dim1.width - w) / 2,
				(dim1.height - h) / 2);		
		about.setVisible(true);
	}

	
	/** Shows the licensing information. */
	public void viewHelpLicense() 
	{
		JFrame licensing = new JFrame("PRIORITY: License Information");
		licensing.setSize(550, 500);
		licensing.setIconImage(Priority.icon);

		String fileName = "/license.html";
		java.net.URL url = this.getClass().getResource(fileName);
	
		if (url != null) 
		{
			JEditorPane licenseArea = null;
			
			try 
			{
			   licenseArea = new JEditorPane(url);
			}
			catch (IOException e) 
			{
				JLabel l = new JLabel("Please visit www.cs.duke.edu for licensing information.");
				l.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
				licensing.getContentPane().add(l);
				licensing.setSize(380, 200);			
			}
			
			licenseArea.setMargin(new Insets(5,5,5,5));
			licenseArea.setCaretPosition(0);
			licenseArea.setContentType("text/html");
			licenseArea.setEditable(false);
			licensing.getContentPane().add(new JScrollPane(licenseArea));

		} 
		else 
		{
			JLabel l = new JLabel("Please visit www.cs.duke.edu for licensing information.");
			l.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
			licensing.getContentPane().add(l);
			licensing.setSize(380, 200);
		}
		
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		licensing.validate();
		licensing.setLocation((dim.width - licensing.getWidth()) / 2, (dim.height - licensing.getHeight()) / 2);
		licensing.setVisible(true);
	}

	
	/** Shows the tips. */
	public void viewHelpTips() 
	{
		JFrame tips = new JFrame("PRIORITY: Tips");
		tips.setSize(550, 380);
		tips.setIconImage(Priority.icon);

		String fileName = "/tips.html";
		java.net.URL url = this.getClass().getResource(fileName);
	
		if (url != null) 
		{
			JEditorPane tipsArea = null;
			
			try 
			{
				tipsArea = new JEditorPane(url);
			}
			catch (IOException e) 
			{
				JLabel l = new JLabel("No tips available.");
				l.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
				tips.getContentPane().add(l);
				tips.setSize(150, 100);			
			}
			
			tipsArea.setMargin(new Insets(5,5,5,5));
			tipsArea.setCaretPosition(0);
			tipsArea.setContentType("text/html");
			tipsArea.setEditable(false);
			tips.getContentPane().add(new JScrollPane(tipsArea));

		} 
		else 
		{
			JLabel l = new JLabel("No tips available.");
			l.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
			tips.getContentPane().add(l);
			tips.setSize(150, 100);
		}
		Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
		tips.validate();
		tips.setLocation((dim.width - tips.getWidth()) / 2, (dim.height - tips.getHeight()) / 2);
		tips.setVisible(true);
	}
}