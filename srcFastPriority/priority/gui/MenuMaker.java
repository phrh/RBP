package priority.gui;
import java.awt.Dimension;
import java.awt.event.KeyEvent;
import java.awt.event.ActionEvent;
import java.io.File;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import javax.swing.SwingConstants;
import javax.swing.JMenuItem;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.BorderFactory;
import javax.swing.AbstractAction;
import javax.swing.ButtonGroup;
/*import javax.swing.JFileChooser;*/
import priority.Strings;

/** MenuMaker - the class creates the main
 * menu for Priority.
 * @author raluca
 * some code borrowed from Matt Edwards, Jason Bosko
 */
class MenuMaker {
	
	private static int menuItemWidth = 50;
	private static int menuItemHeight = 20;
	public static JMenuItem helpMenuLicense, helpMenuReadme, helpMenuAbout, helpMenuTips;
	public static JMenuItem viewMenuViewPSSM, viewMenuViewPositions;
	public static JRadioButtonMenuItem priorsMenuSinglePrior, priorsMenuMultiplePriors;
	public static ButtonGroup group;
	
	
	public static void makeMenus(final JFrame parent) 
	{
		JMenuBar bar = new JMenuBar();
		/*bar.add(makeFileMenu(parent));*/
		bar.add(makePriorsMenu(parent));
		bar.add(makeViewMenu(parent));
		bar.add(makeHelpMenu(parent));
		parent.setJMenuBar(bar);
	}
	
	private static JMenu makePriorsMenu(final JFrame parent) 
	{
		JMenu priorsMenu = new JMenu(Strings.getString("priorsMenu"));
		priorsMenu.setPreferredSize(new Dimension(menuItemWidth, menuItemHeight));
		priorsMenu.setHorizontalAlignment(SwingConstants.CENTER);
		priorsMenu.setMnemonic(KeyEvent.VK_P);
		
		priorsMenuSinglePrior = new JRadioButtonMenuItem();
		priorsMenuSinglePrior.setBorder(BorderFactory.createEmptyBorder(3,3,3,3));
		priorsMenuSinglePrior.setMnemonic(KeyEvent.VK_S);
		priorsMenu.add(priorsMenuSinglePrior);
		priorsMenuSinglePrior.setAction(new AbstractAction(Strings.getString("priorsMenuSinglePrior")) {
			private static final long serialVersionUID = 1;
			public void actionPerformed(ActionEvent evt) {					    
				((MainWindow)parent).switchToSinglePrior();
			}
		});
		
		priorsMenuMultiplePriors = new JRadioButtonMenuItem();
		priorsMenuMultiplePriors.setBorder(BorderFactory.createEmptyBorder(3,3,3,3));
		priorsMenuMultiplePriors.setMnemonic(KeyEvent.VK_M);
		priorsMenu.add(priorsMenuMultiplePriors);
		priorsMenuMultiplePriors.setAction(new AbstractAction(Strings.getString("priorsMenuMultiplePriors")) {
			private static final long serialVersionUID = 1;
			public void actionPerformed(ActionEvent evt) {					    
				((MainWindow)parent).switchToMultiplePriors();
			}
		});
		
		group = new ButtonGroup();
		group.add(priorsMenuSinglePrior);
		group.add(priorsMenuMultiplePriors);
		
		return priorsMenu;
	}
	
	
	private static JMenu makeViewMenu(final JFrame parent) 
	{
		JMenu viewMenu = new JMenu(Strings.getString("viewMenu"));
		viewMenu.setPreferredSize(new Dimension(menuItemWidth, menuItemHeight));
		viewMenu.setHorizontalAlignment(SwingConstants.CENTER);
		viewMenu.setMnemonic(KeyEvent.VK_V);
		
		viewMenuViewPSSM = new JMenuItem();
		viewMenuViewPSSM.setBorder(BorderFactory.createEmptyBorder(3,3,3,3));
		viewMenuViewPSSM.setMnemonic(KeyEvent.VK_S);
		viewMenu.add(viewMenuViewPSSM);
		viewMenuViewPSSM.setAction(new AbstractAction(Strings.getString("viewMenuViewPSSM")) {
			private static final long serialVersionUID = 1;
			public void actionPerformed(ActionEvent evt) {					    
				((MainWindow)parent).viewFinalPSSM();
			}
		});
		
		viewMenuViewPositions = new JMenuItem();
		viewMenuViewPositions.setBorder(BorderFactory.createEmptyBorder(3,3,3,3));
		viewMenuViewPositions.setMnemonic(KeyEvent.VK_O);
		viewMenu.add(viewMenuViewPositions);
		viewMenuViewPositions.setAction(new AbstractAction(Strings.getString("viewMenuViewPositions")) {
			private static final long serialVersionUID = 1;
			public void actionPerformed(ActionEvent evt) {					    
				((MainWindow)parent).viewFinalMotifPositions();
			}
		});
		
		return viewMenu;
	}
	
	
	private static JMenu makeHelpMenu(final JFrame parent) 
	{
		JMenu helpMenu = new JMenu(Strings.getString("helpMenu"));
		helpMenu.setPreferredSize(new Dimension(menuItemWidth, menuItemHeight));
		helpMenu.setHorizontalAlignment(SwingConstants.CENTER);			
		helpMenu.setMnemonic('H');
		
		helpMenuReadme = new JMenuItem();
		helpMenuReadme.setBorder(BorderFactory.createEmptyBorder(3,3,3,3));
		helpMenuReadme.setMnemonic(KeyEvent.VK_M);
		helpMenuReadme.setAction(new AbstractAction(Strings.getString("helpMenuReadme")) {
			private static final long serialVersionUID = 1;
			public void actionPerformed(ActionEvent evt) {					    
				((MainWindow)parent).viewHelpReadme();
			}
		});	
		
		helpMenuTips = new JMenuItem();
		helpMenuTips.setBorder(BorderFactory.createEmptyBorder(3,3,3,3));
		helpMenuTips.setMnemonic(KeyEvent.VK_T);
		helpMenuTips.setAction(new AbstractAction(Strings.getString("helpMenuTips")) {
			private static final long serialVersionUID = 1;
			public void actionPerformed(ActionEvent evt) {					    
				((MainWindow)parent).viewHelpTips();
			}
		});	

		helpMenuLicense = new JMenuItem();
		helpMenuLicense.setBorder(BorderFactory.createEmptyBorder(3,3,3,3));
		helpMenuLicense.setMnemonic(KeyEvent.VK_L);
		helpMenuLicense.setAction(new AbstractAction(Strings.getString("helpMenuLicense")) {
			private static final long serialVersionUID = 1;
			public void actionPerformed(ActionEvent evt) {					    
				((MainWindow)parent).viewHelpLicense();
			}
		});		


		helpMenuAbout = new JMenuItem();
		helpMenuAbout.setBorder(BorderFactory.createEmptyBorder(3,3,3,3));
		helpMenuAbout.setMnemonic(KeyEvent.VK_A);
		helpMenuAbout.setAction(new AbstractAction(Strings.getString("helpMenuAbout")) {
			private static final long serialVersionUID = 1;
			public void actionPerformed(ActionEvent evt) {					    
				((MainWindow)parent).viewHelpAbout();
			}
		});
		helpMenu.add(helpMenuLicense);
		helpMenu.add(helpMenuReadme);
		helpMenu.add(helpMenuTips);
		helpMenu.addSeparator();
		helpMenu.add(helpMenuAbout);
		return helpMenu;
	}
}


class ParamFileFilter extends javax.swing.filechooser.FileFilter 
{
    
	public boolean accept(File file) 
	{
        return (file.getName()).endsWith("." + Strings.getString("paramFilesExtension"));
    }
	
    public String getDescription() 
    {
        return "*." + Strings.getString("paramFilesExtension");
    }
}