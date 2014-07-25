package priority.gui;
import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.Dimension;
import javax.swing.JRadioButton;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.JSpinner;
import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.SpinnerNumberModel;
import javax.swing.SpinnerModel;
import javax.swing.DefaultListModel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import priority.Parameters;
import priority.Strings;


/**
 * ParamPanel - a panel for setting the parameters. 
 * @author raluca
 */
class ParamPanel extends JPanel implements ActionListener 
{
	static final long serialVersionUID = 1;
	
	JPanel north, center, pathPanel, numberPanel;
	
	/* in north panel */
	JTextField tfName;
	FilePanel tfFilePanel;
	JRadioButton fileTF, individualTF;
	ButtonGroup radioButtonGroup;
	
	/* in center->pathPanel */
	JLabel pathData, pathBkgr, pathOutput;
	FilePanel dataFilePanel, bkgrFilePanel, outputFilePanel;

	/* in center->center */
	PriorFilesPanel priorFilesPanel;
	SinglePriorFilePanel singlePriorFilePanel;
	
	/* in center->numberPanel */
	JSpinner wsizeSpinner, iterationsSpinner, trialsSpinner, bkgrOrderSpinner;
	
	
	/** Creates and initilizes the components. */
	public void init() 
	{
		this.setLayout(new BorderLayout());
		this.setBorder(BorderFactory.createEmptyBorder(5,0,0,5));
		
		/* the north panel */
		north = new JPanel(new GridLayout(4, 1, 0, 5));
		north.setBorder(BorderFactory.createCompoundBorder( 
				BorderFactory.createCompoundBorder(
					BorderFactory.createEmptyBorder(1,1,0,1),
					BorderFactory.createLoweredBevelBorder()),
				BorderFactory.createCompoundBorder(
						BorderFactory.createRaisedBevelBorder(),
						BorderFactory.createEmptyBorder(0,10,10,10)))
				);
		north.setPreferredSize(new Dimension(190,125));
		
		fileTF = new JRadioButton(Strings.getString("fileTF"));
		fileTF.setPreferredSize(new Dimension(60,20));
		fileTF.setBorder(BorderFactory.createEmptyBorder(8,0,3,0));		
		fileTF.addActionListener(this);
		
		individualTF = new JRadioButton(Strings.getString("individualTF"));
		individualTF.setPreferredSize(new Dimension(60,20));
		individualTF.setBorder(BorderFactory.createEmptyBorder(8,0,3,0));
		individualTF.addActionListener(this);
		radioButtonGroup = new ButtonGroup();
		radioButtonGroup.add(fileTF);
		radioButtonGroup.add(individualTF);
		
		tfFilePanel = new FilePanel(Parameters.tf_path, false);
		tfName = new JTextField();
		tfName.setPreferredSize(new Dimension(60,20));
		tfName.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createLoweredBevelBorder(),
				BorderFactory.createEmptyBorder(0,5,0,5)));
		tfName.setText(Parameters.tf_name);
		
		if (Parameters.individualTF == false) {
			fileTF.setSelected(true);
			tfName.setEnabled(false);
		} 
		else {
			individualTF.setSelected(true);
			tfFilePanel.setEnabled(false);
		}
		
		north.add(individualTF);
		north.add(tfName);
		north.add(fileTF);
		north.add(tfFilePanel);		
		this.add(north, BorderLayout.NORTH);
		
		/* the center panel */
		center = new JPanel(new BorderLayout());
		this.add(center, BorderLayout.CENTER);
		
		/* the center->north panel (pathPanel) */
		pathPanel = new JPanel(new GridLayout(6, 1, 0, 0));	
		pathPanel.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createCompoundBorder(
					BorderFactory.createEmptyBorder(7,1,1,1),
					BorderFactory.createLoweredBevelBorder()),
				BorderFactory.createCompoundBorder(
						BorderFactory.createRaisedBevelBorder(),
						BorderFactory.createEmptyBorder(0,10,5,10)))
				);
		pathPanel.setPreferredSize(new Dimension(70,170));
			
		dataFilePanel = new FilePanel(Parameters.fname_path, true);
		bkgrFilePanel = new FilePanel(Parameters.back_file, false);
		outputFilePanel = new FilePanel(Parameters.path_output, true);

		pathData = new JLabel(Strings.getString("pathData"));
		pathData.setBorder(BorderFactory.createEmptyBorder(15,0,10,0));
		pathBkgr = new JLabel(Strings.getString("pathBkgr"));
		pathBkgr.setBorder(BorderFactory.createEmptyBorder(15,0,10,0));
		pathOutput = new JLabel(Strings.getString("pathOutput"));
		pathOutput.setBorder(BorderFactory.createEmptyBorder(15,0,10,0));

		pathPanel.add(pathData);
		pathPanel.add(dataFilePanel);
		pathPanel.add(pathBkgr);
		pathPanel.add(bkgrFilePanel);
		pathPanel.add(pathOutput);
		pathPanel.add(outputFilePanel);
		center.add(pathPanel, BorderLayout.NORTH);		
		
		
		/* the center->center panel */
		priorFilesPanel = new PriorFilesPanel();
		priorFilesPanel.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createCompoundBorder(
					BorderFactory.createEmptyBorder(7,1,1,1),
					BorderFactory.createLoweredBevelBorder()),
				BorderFactory.createCompoundBorder(
						BorderFactory.createRaisedBevelBorder(),
						BorderFactory.createEmptyBorder(10,10,5,10)))
				);
		
		singlePriorFilePanel = new SinglePriorFilePanel();
		singlePriorFilePanel.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createCompoundBorder(
					BorderFactory.createEmptyBorder(7,1,1,1),
					BorderFactory.createLoweredBevelBorder()),
				BorderFactory.createCompoundBorder(
						BorderFactory.createRaisedBevelBorder(),
						BorderFactory.createEmptyBorder(5,5,5,5)))
				);
		
		if (Parameters.multiple_priors == false)
			center.add(singlePriorFilePanel, BorderLayout.CENTER);
		else 
			center.add(priorFilesPanel, BorderLayout.CENTER);
		
		
		/* the center->south panel (numberPanel) */
		numberPanel = new JPanel(new BorderLayout());
		JPanel numberPanelCenter = new JPanel(new GridLayout(4, 1, 3, 5));
		numberPanel.add(numberPanelCenter, BorderLayout.CENTER);
		JPanel numberPanelEast = new JPanel(new GridLayout(4, 1, 3, 5));
		numberPanelEast.setBorder(BorderFactory.createEmptyBorder(0,5,0,5));
		
		numberPanel.add(numberPanelEast, BorderLayout.EAST);
		numberPanel.setBorder(BorderFactory.createCompoundBorder(
				BorderFactory.createCompoundBorder(
					BorderFactory.createEmptyBorder(7,1,1,1),
					BorderFactory.createLoweredBevelBorder()),
				BorderFactory.createCompoundBorder(
						BorderFactory.createRaisedBevelBorder(),
						BorderFactory.createEmptyBorder(5,0,5,0)))
				);
		
		wsizeSpinner = ParamPanel.createAndAddSpinner(numberPanelCenter, numberPanelEast,
				Strings.getString("wsizeSpinner"), Parameters.wsize, 
				Parameters.wsizeMin, Parameters.wsizeMax);
		bkgrOrderSpinner = ParamPanel.createAndAddSpinner(numberPanelCenter, numberPanelEast,
				Strings.getString("bkgrOrderSpinner"), Parameters.bkgrOrder, 
				Parameters.bkgrOrderMin, Parameters.bkgrOrderMax);

		trialsSpinner = ParamPanel.createAndAddSpinner(numberPanelCenter, numberPanelEast,
				Strings.getString("trialsSpinner"), Parameters.trials, 1, Integer.MAX_VALUE);
		iterationsSpinner = ParamPanel.createAndAddSpinner(numberPanelCenter, numberPanelEast, 
				Strings.getString("iterationsSpinner"), Parameters.iter, 1, Integer.MAX_VALUE);
		
		center.add(numberPanel, BorderLayout.SOUTH);
	}
	
	
	/** Creates spinner, sets spinner properties and adds it to the panel. */
	private static JSpinner createAndAddSpinner(JPanel panelLabel, JPanel panel, 
			String label, int value, int min, int max) {
		SpinnerNumberModel model = new SpinnerNumberModel(value, min, max, 1);
		
		JSpinner spinner = new JSpinner(model);
		spinner.setPreferredSize(new Dimension(80,22));
		spinner.setBorder(BorderFactory.createLoweredBevelBorder());
		
		JLabel ll = new JLabel(label,JLabel.RIGHT);
		ll.setLabelFor(spinner);
		((JSpinner.DefaultEditor)(spinner.getEditor())).getTextField().
			setBorder(BorderFactory.createEmptyBorder(0,5,0,5));
	
		panelLabel.add(ll);
		panel.add(spinner);
		return spinner;
	}

	
	/** Implements the actionPerformed function (from ActionListener) */
	public void actionPerformed(ActionEvent ev) {
		if (ev.getSource() == fileTF) {
			tfName.setEnabled(false);
			tfFilePanel.setEnabled(true);
 			return;
 		}
		if (ev.getSource() == individualTF) {
			tfName.setEnabled(true);
			tfFilePanel.setEnabled(false);
 			return;
 		}
	}	
	

	public void switchToSinglePrior() {
		center.remove(priorFilesPanel);
		center.add(singlePriorFilePanel, BorderLayout.CENTER);   
		Parameters.multiple_priors = false;
	}
	
	public void switchToMultiplePriors() {
		center.remove(singlePriorFilePanel);
		center.add(priorFilesPanel, BorderLayout.CENTER);
		Parameters.multiple_priors = true;
	}
	
	
	/** Gets the values for all parameters from the widgets 
	 *  in this panel, checks them and copies them into the 
	 *  static members of the class Parameters. */ 
	public String setParameterValues() {
		
		/* the number of iterations and number of trials are valid for sure */
		/* we also set wsize and bkgrOrder (whether or not they are valid will be checked later) */
		try {
			Parameters.trials = Integer.parseInt((String)
					((JSpinner.DefaultEditor)trialsSpinner.getEditor()).getTextField().getText());
			Parameters.iter = ((Number)((SpinnerModel)iterationsSpinner.getModel()).getValue()).intValue();
			Parameters.wsize = Integer.parseInt((String)
					((JSpinner.DefaultEditor)wsizeSpinner.getEditor()).getTextField().getText());
			Parameters.bkgrOrder = Integer.parseInt((String)
					((JSpinner.DefaultEditor)bkgrOrderSpinner.getEditor()).getTextField().getText());
		}
		catch(Exception e) {
			System.out.println(e);/* it shouldn't get here (ever) */
		}

		
	    /* *********************************************************************************** */
		/* the path for the output files: must be a writable directory */
		String path = outputFilePanel.getText();
		File dir = new File(path);
	    if (!dir.exists()) 
	    	return "Error: the directory for the output files (\"" + path + "\") does not exist!";
		if (!dir.isDirectory()) 
	    	return "Error: the path for the output files (\"" + path + "\") is not a valid directory!";
	    if (!dir.canWrite()) 
	    	return "Error: the directory for the output files (\"" + path + "\") is not writable!";
		Parameters.path_output = path;
		
		/* the bkgr model file: must be a readable file + the content will be checked later */
		path = bkgrFilePanel.getText();
		File file = new File(path);
	    if (!file.exists()) 
	    	return "Error: the background model file (\"" + path + "\") does not exist!";
	    if (!file.isFile()) 
	    	return "Error: the background model file (\"" + path + "\") is not a valid file!";
	    if (!file.canRead()) 
	    	return "Error: the background model file (\"" + path + "\") is not readable!";
		Parameters.back_file = path;
		
		/* the path for the FASTA data files: must be a readable directory + check content later */
		path = dataFilePanel.getText();
		dir = new File(path);
	    if (!dir.exists()) 
	    	return "Error: the directory for the FASTA files (\"" + path + "\") does not exist!";
		if (!dir.isDirectory()) 
	    	return "Error: the path for the FASTA files (\"" + path + "\") is not a valid directory!";
	    if (!dir.canRead()) 
	    	return "Error: the directory for the FASTA files (\"" + path + "\") is not readable!";
		Parameters.fname_path = path;

		
		/* *********************************************************************************** */
		/* next we check the content of the bkgr model file: it must have 4 or 4+4^2 or
		 * 4+4^2+4^3 or ...+4^k lines, each line must be between 0 and 1 and each group of 
		 * four lines must add to 1 (or very close to 1) + the order must be <(k-1) */ 
		/* the bkgr file name and the order are already set in the class Parameters,
		 * so we don't need to send them !!!! */
		String err = Parameters.readBackground();
		if (err.compareTo("") != 0)
			return err;
		System.out.println("Background model... OK");
		

		/* setting the TF names */
		if (individualTF.isSelected()) { /*we apply the alg for a TF only */
			Parameters.tf_name = (tfName.getText()).trim();
			if (Parameters.tf_name.length() < 1) {
				tfName.requestFocus();
				return "Error: the TF name field is empty!";
			}
			Parameters.individualTF = true;
			Parameters.tf_names = null;
			Parameters.tf_names = new String[1];
			Parameters.tf_names[0] = Parameters.tf_name;
		}
		else { /* we apply the alg to all the TFs in a file */
			/* the TF names file must be a readable file */
			path = tfFilePanel.getText();
			file = new File(path);
		    if (!file.exists()) 
		    	return "Error: the file with the TF names (\"" + path + "\") does not exist!";
		    if (!file.isFile()) 
		    	return "Error: the file with the TF names (\"" + path + "\") is not a valid file!";
		    if (!file.canRead()) 
		    	return "Error: the file with the TF names (\"" + path + "\") is not readable!";
			Parameters.tf_path = path;
			Parameters.individualTF = false;
			err = Parameters.read_TFnames();
			if (err.compareTo("") != 0)
				return err;
		}
		
		/* next we have to check that for every TF there is a file in the data path
		 * with the same name as the TF */
		for (int i=0; i<Parameters.tf_names.length; i++) {
			path = Parameters.fname_path + "/" + Parameters.tf_names[i] + ".fasta";
			file = new File(path);
		    if (!file.exists()) 
		    	return "Error: the data file \"" + path + "\"\ncorresponding to the TF \"" +
		    	Parameters.tf_names[i] + "\" does not exist!";
		    if (!file.isFile()) 
		    	return "Error: the data file \"" + path + "\"\ncorresponding to the TF \"" +
		    	Parameters.tf_names[i] + "\" is not a valid file!";
		    if (!file.canRead()) 
		    	return "Error: the data file \"" + path + "\"\ncorresponding to the TF \"" +
		    	Parameters.tf_names[i] + "\" is not a readable file!";
		}
		
		
		/* *********************************************************************************** */
		/* SET THE PRIOR INFORMATION */
		
		/* first check single/multiple files */
		if (Parameters.multiple_priors == false) { /* uniform or single prior */
			Parameters.putative_class = -1;
			if (this.singlePriorFilePanel.uniformPrior.isSelected()) {
				Parameters.otherclass = true;
				Parameters.prior_dirs = new String[0];
			}
			else {
				Parameters.otherclass = false;
				Parameters.prior_dirs = new String[1];
				Parameters.prior_dirs[0] = this.singlePriorFilePanel.priorFilePanel.getText();
			}
		}
		else { /* multiple priors */
			/* set the putative class and the "other class"-flag */
			if (priorFilesPanel.putativeClassCheckbox.isSelected()) 
				Parameters.putative_class = ((DefaultListModel)priorFilesPanel.classList.getModel()).
					indexOf(priorFilesPanel.putativeClassText.getText()); 
			else
				Parameters.putative_class = -1;
			Parameters.otherclass = priorFilesPanel.otherClassCheckbox.isSelected(); 

	    	/* now the class names will be copied in Parameters.prior_dirs */
		    int n;
		    if (Parameters.otherclass)
		    	n = priorFilesPanel.classListModel.size() - 1;
		    else 
		    	n = priorFilesPanel.classListModel.size();
		    Parameters.prior_dirs = new String[n];
		    for (int i=0;i<n;i++)
		    	Parameters.prior_dirs[i] = (String)priorFilesPanel.classListModel.getElementAt(i);		    
	    }
		
	    /* Each entry in Parameters.prior_dirs (if such entries exist) must be the name 
	     * of a readable directory. */	
		for (int i=0; i<Parameters.prior_dirs.length; i++) {
			dir = new File(Parameters.prior_dirs[i]);
			if ((!dir.exists()) || (!dir.isDirectory()) || (!dir.canRead()))
		    	return "Error: the prior directory: \"" + Parameters.prior_dirs[i] + 
		    		"\" does not exist or is not readable!";			
		}

		if (err.compareTo("") != 0)
			return err;
		System.out.println("Prior files (" + Parameters.prior_dirs.length + ")... OK");

		return "";
	}
}