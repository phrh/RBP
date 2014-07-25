package priority.gui;
import java.awt.BorderLayout;
import java.awt.Dimension;
import javax.swing.JTextField;
import javax.swing.JPanel;
import javax.swing.JButton;
import javax.swing.BorderFactory;
import javax.swing.JFileChooser;
import java.awt.event.ActionEvent; 
import java.awt.event.ActionListener;

/**
 * FilePanel - a panel that incorporates a 
 * file chooser and a text field. 
 * @author raluca
 */
class FilePanel extends JPanel implements ActionListener
{
	static final long serialVersionUID = 1;
	private JTextField fileNameField;
	private JButton choose;
	private boolean dir;
	private String currentParentName = ".";
	
	/** Constructor */
	public FilePanel(String fileName, boolean directory) {
		dir = directory;
		init(fileName);
		java.io.File file = new java.io.File(fileName);
		if (file.exists()) {
			if ((dir == true) && file.isDirectory())
				currentParentName = file.getParent();
			if ((dir == false) && file.isFile())
				currentParentName = file.getParent();
		}
		file = null;
	}
	
	/** Creates and initilizes the components. */
	private void init(String fileName) {
		BorderLayout borderLayout = new BorderLayout();
		borderLayout.setHgap(10);
		this.setLayout(borderLayout);

		choose = new JButton("...");
		String os = "unknown";
		try { os = System.getProperty("os.name");
		} catch (SecurityException e) {}
		if(os.indexOf("Mac") >= 0) 
			choose.setBorder(BorderFactory.createRaisedBevelBorder());			
		choose.setPreferredSize(new Dimension(30,20));
		choose.addActionListener(this);	
		this.add(choose, BorderLayout.EAST);
		
		fileNameField = new JTextField();
		fileNameField.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createLoweredBevelBorder(),
				BorderFactory.createEmptyBorder(0,5,0,5)));
		fileNameField.setText(fileName);
		fileNameField.setPreferredSize(new Dimension(60,15));
		fileNameField.addActionListener(this);
		this.add(fileNameField, BorderLayout.CENTER);
	}
	
	/** Returns the file name. */
	public String getText()
	{
		return fileNameField.getText();
	}
	
    /** Handles the button click. */
	public void actionPerformed(ActionEvent ev) {
		if (ev.getSource() == choose) {
			JFileChooser ch = new JFileChooser(this.currentParentName);
			if (dir == true) 
				ch.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			int r = ch.showOpenDialog(this);
			if (r == JFileChooser.APPROVE_OPTION) {
			    fileNameField.setText(ch.getSelectedFile().getAbsolutePath() );
			}
		}
	}
	    
    /** Enables or disables the panel. */
    public void setEnabled(boolean state) {
   		fileNameField.setEnabled(state);
   		choose.setEnabled(state);
    }
}