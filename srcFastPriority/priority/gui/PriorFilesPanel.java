package priority.gui;
import java.awt.*;
import javax.swing.JPanel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JTextField;
import javax.swing.DefaultListModel;
import javax.swing.JList;
import javax.swing.JFileChooser;
import javax.swing.BorderFactory;
import javax.swing.JScrollPane;
import javax.swing.JOptionPane;
import javax.swing.ListSelectionModel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemListener;
import java.awt.event.ItemEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseEvent;
import java.awt.event.KeyListener;
import java.awt.event.KeyEvent;
import java.io.File;
import priority.Parameters;
import priority.Strings;
import java.util.ArrayList;


/**
 * PriorFilesPanel - a panel that incorporates a 
 * directory chooser, a list, a combobox and 
 * and 2 checkboxes, used to set the class prior dirs. 
 * @author raluca
 */
class PriorFilesPanel extends JPanel implements ActionListener, ItemListener, MouseListener, KeyListener
{
	static final long serialVersionUID = 1;
	public String currentDir;
	
	private JButton priorDirsButton;
	private JPanel center, smallPanel; 
	JCheckBox otherClassCheckbox, putativeClassCheckbox;
	JTextField putativeClassText;
	DefaultListModel classListModel;
	JList classList;
	String other_class = "other prior-type";
	String no_putative_class ="no putative prior-type";
	
	
	/** Constructor */
	public PriorFilesPanel() 
	{
		init();
	}
	
	
	/** Creates and initilizes the components. */
	private void init() 
	{
		BorderLayout borderLayout = new BorderLayout();  
		borderLayout.setVgap(10);
		setLayout(borderLayout);
		setPreferredSize(new Dimension(100, 220));

		if (Parameters.prior_dirs.length == 0)
		{
			currentDir = ".";
		}
		else
		{
			currentDir = Parameters.prior_dirs[0];
		}
		
		
		priorDirsButton = new JButton(Strings.getString("priorFilesButton"));
		String os = "unknown";
		
		try 
		{ 
			os = System.getProperty("os.name");
		} 
		catch (SecurityException e) {}
		
		if(os.indexOf("Mac") >= 0) 
		{
			priorDirsButton.setBorder(BorderFactory.createRaisedBevelBorder());
		}
		
		
		priorDirsButton.setPreferredSize(new Dimension(100,25));
		priorDirsButton.addActionListener(this);	
		this.add(priorDirsButton, BorderLayout.NORTH);
				
		center = new JPanel(new BorderLayout());
		this.add(center, BorderLayout.CENTER);
		
		classListModel = new DefaultListModel();
		for (int i=0; i<Parameters.prior_dirs.length; i++)
		{
			classListModel.addElement(Parameters.prior_dirs[i]);
		}
		
		if (Parameters.otherclass)
		{
			classListModel.addElement(other_class);
		}
		
		classList = new JList(classListModel);
		classList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		classList.setLayoutOrientation(JList.VERTICAL);
		classList.setVisibleRowCount(0);
		
		if (classListModel.size() > 0)
		{
			classList.setSelectedIndex(0);
		}
		
		JScrollPane listScroller = new JScrollPane(classList, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		
		listScroller.setBorder(BorderFactory.createLoweredBevelBorder());
		classList.addMouseListener(this);
		classList.addKeyListener(this);
		center.add(listScroller, BorderLayout.CENTER);
		
		smallPanel = new JPanel(new GridLayout(2,1,0,0));		
		putativeClassCheckbox = new JCheckBox(Strings.getString("putativeClassCheckbox"));
		putativeClassCheckbox.setBorder(BorderFactory.createEmptyBorder(11,5,1,5));
		putativeClassCheckbox.addItemListener(this);

		putativeClassText = new JTextField();
		putativeClassText.setPreferredSize(new Dimension(60,15));
		putativeClassText.setBorder(
				BorderFactory.createCompoundBorder(
						BorderFactory.createEmptyBorder(3,5,0,5),
						BorderFactory.createCompoundBorder(
								BorderFactory.createLoweredBevelBorder(),
								BorderFactory.createEmptyBorder(0,5,0,5))));
		putativeClassText.setEditable(false);
		
		if ((Parameters.putative_class >= 0) && (Parameters.putative_class < classListModel.size())) 
		{
			putativeClassText.setText((String)classListModel.elementAt(Parameters.putative_class));
			putativeClassText.setEnabled(true);
			classList.setSelectedIndex(Parameters.putative_class);
			putativeClassCheckbox.setSelected(true);
			putativeClassText.setBackground(smallPanel.getBackground());
		}
		else 
		{
			putativeClassText.setText(no_putative_class);
			putativeClassText.setEnabled(false);
			putativeClassCheckbox.setSelected(false);
			putativeClassText.setBackground(smallPanel.getBackground());
		}
		
		smallPanel.add(putativeClassCheckbox);
		smallPanel.add(putativeClassText);
		center.add(smallPanel, BorderLayout.SOUTH);
		
		otherClassCheckbox = new JCheckBox(Strings.getString("otherClassCheckbox"));
		otherClassCheckbox.setSelected(Parameters.otherclass);
		otherClassCheckbox.setBorder(BorderFactory.createEmptyBorder(0,5,4,5));
		otherClassCheckbox.setPreferredSize(new Dimension(100, 20));
		otherClassCheckbox.addItemListener(this);
		this.add(otherClassCheckbox, BorderLayout.SOUTH);
	}

	
    /** Handles the button click. */
	public void actionPerformed(ActionEvent ev) 
	{
		if (ev.getSource() == this.priorDirsButton) 
		{
			JFileChooser ch = new JFileChooser(currentDir);
			ch.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			ch.setMultiSelectionEnabled(true);
			int r = ch.showOpenDialog(this);
			
			if (r == JFileChooser.APPROVE_OPTION) 
			{
				File dirs[] = ch.getSelectedFiles();
				String str; File dir;
				
				/* first take all the valid dir names and put them in an array of strings */
				ArrayList<String> dirNames = new ArrayList<String>(); 
				for (int i=0; i<dirs.length; i++) 
				{
					str = dirs[i].getAbsolutePath();
					dir = new File(str);
				    
					if ((!dir.exists()) || (!dir.isDirectory()) || (!dir.canWrite())) 
					{
				    	String mess =  "Error: " + str + "does not exists or it is not a writable directory!";
				    	System.out.println(mess+"\n");
						JOptionPane.showMessageDialog(null,	str, 
								"Error", JOptionPane.ERROR_MESSAGE);									
				    }
				    else 
				    {
				    	dirNames.add( str );
				    }
				}
				
				/* if the array is not empty we update classList */
				if (dirNames.size() == 0) 
				{
					JOptionPane.showMessageDialog(null,	"No valid directory was selected!", 
							"Error", JOptionPane.ERROR_MESSAGE);									
					return;
				}
					
				/* add the dir names in the list and also in Parameters.prior_dirs */		
				/* the new list will contain the existing dirs plus the newly selected dirs */
				
				String[] temp = new String[Parameters.prior_dirs.length];
				
				for (int i=0; i<Parameters.prior_dirs.length; i++) 
				{
					temp[i] = "" + Parameters.prior_dirs[i];
				}
				
				Parameters.prior_dirs = null;
				Parameters.prior_dirs = new String[temp.length + dirNames.size()];
				
				for (int i=0; i<temp.length; i++) 
				{
					Parameters.prior_dirs[i] = "" + temp[i];
				}
				
				for (int i=0; i<dirNames.size(); i++) 
				{
					Parameters.prior_dirs[temp.length+i] = (String)dirNames.get(i);
					classListModel.addElement(Parameters.prior_dirs[temp.length+i]);
				}
				
				currentDir = dirs[0].getParentFile().getAbsolutePath();
			}
		}
	}
	
    
	/** Implements the itemStateChanged function (from ItemListener) */
	public void itemStateChanged (ItemEvent ev) 
	{
		if (ev.getSource() == putativeClassCheckbox) 
		{
			if (putativeClassCheckbox.isSelected()) 
			{
				/* the selected item in the classList becomes the putative class */
				int index = classList.getSelectedIndex();
				
				if (index < 0) 
				{
					/* no item is selected. We select the first one, if it exists. */
					if (classListModel.size() > 0) 
					{
						classList.setSelectedIndex(0);
						JOptionPane.showMessageDialog(null, "No prior-type was selected, so the first prior-type\n"+
							"in the list became the putative prior-type.", "Warning", JOptionPane.WARNING_MESSAGE);
						putativeClassText.setText((String)classList.getModel().getElementAt(0));
						putativeClassText.setEnabled(true);
						putativeClassText.setBackground(smallPanel.getBackground());
					}
					else 
					{
						/* it should not get here */
						JOptionPane.showMessageDialog(null, "There are no prior-types in the list!", 
								"Error", JOptionPane.ERROR_MESSAGE); 
						putativeClassCheckbox.setSelected(false);
						putativeClassText.setText(no_putative_class);
						putativeClassText.setEnabled(false);
						putativeClassText.setBackground(smallPanel.getBackground());
					}
				}
				else 
				{
					putativeClassText.setText((String)classList.getModel().getElementAt(
							classList.getSelectedIndex()));
					putativeClassText.setEnabled(true);
					putativeClassText.setBackground(smallPanel.getBackground());					
				}
			}
			else 
			{
				putativeClassText.setText(no_putative_class);
				putativeClassText.setEnabled(false);
				putativeClassText.setBackground(smallPanel.getBackground());
			}
			
			return;
		}
		
		if (ev.getSource() == otherClassCheckbox) 
		{
			if (otherClassCheckbox.isSelected())
			{
				classListModel.addElement(other_class);
			}
			else 
			{ 
				if (classListModel.size() == 1)
				{
					classListModel.removeElement(other_class);
					otherClassCheckbox.setSelected(true);
					classList.setSelectedIndex(0);
					return;
				}
				
				int index = classListModel.indexOf(other_class);
				
				if ((classList.getSelectedIndex() == index) /* other_class was selected */
				    && putativeClassCheckbox.isSelected())
				{
					putativeClassCheckbox.setSelected(false);
				}
				
				classListModel.removeElement(other_class);
			}
			return;
		}
	}

	
	/** Handles the mouse events generated by the list of class names. */
	public void mouseClicked(MouseEvent ev) 
	{
		if (ev.getClickCount() == 2) 
		{
			putativeClassText.setText((String)classList.getModel().getElementAt(
					classList.getSelectedIndex()));
			putativeClassText.setEnabled(true);
			putativeClassText.setBackground(smallPanel.getBackground());
			putativeClassCheckbox.setSelected(true);
		}
	}
	
	
	public void mousePressed(MouseEvent ev) {}
	public void mouseReleased(MouseEvent ev) {}	
	public void mouseEntered(MouseEvent ev) {}
	public void mouseExited(MouseEvent ev) {}

	public void keyTyped(KeyEvent ev) {}
	public void keyPressed(KeyEvent ev) {}
	
	public void keyReleased(KeyEvent ev) 
	{
		if ((ev.getKeyCode() == KeyEvent.VK_DELETE) || (ev.getKeyCode() == KeyEvent.VK_BACK_SPACE)) 
		{
			if (classListModel.size() == 1) 
			{
				JOptionPane.showMessageDialog(this, "Error: The prior-types list cannot remain empty!",
						"Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			int index = classList.getSelectedIndex();
			
			if (classList.getModel().getElementAt(index) == other_class) 
			{
				otherClassCheckbox.setSelected(false);
				classList.setSelectedIndex(0);
				return;
			}
			
			if (putativeClassCheckbox.isSelected() && 
			   putativeClassText.getText().compareTo((String)classList.getModel().getElementAt(index)) == 0)
		    { /* the putative class must be deleted */
				putativeClassCheckbox.setSelected(false);
				classListModel.remove(index);
		    }
			else
			{
				classListModel.remove(index);
			}
			
			classList.setSelectedIndex(0);
		}
	}
}