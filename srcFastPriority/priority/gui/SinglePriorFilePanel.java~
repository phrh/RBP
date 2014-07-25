package priority.gui;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;
import javax.swing.ButtonGroup;
import javax.swing.JRadioButton;
import javax.swing.JPanel;
import javax.swing.BorderFactory;
import priority.Parameters;
import priority.Strings;
import java.awt.event.ActionEvent; 
import java.awt.event.ActionListener;


/**
 * SinglePriorFilePanel - a panel that allows the user  
 * to choose between uniform prior and a positional prior
 * @author raluca
 */
public class SinglePriorFilePanel extends JPanel implements ActionListener
{
	static final long serialVersionUID = 1;
	
	FilePanel priorFilePanel;
	JRadioButton uniformPrior, singlePrior;
	ButtonGroup radioButtonGroup;

	/** Constructor */
	public SinglePriorFilePanel() {
		init();
	}

	
	/** Creates and initilizes the components. */
	public void init() 
	{
		this.setLayout(new BorderLayout());

		/* the panel */
		JPanel panel = new JPanel(new GridLayout(3, 1, 0, 3));
		panel.setBorder(BorderFactory.createEmptyBorder(0,5,0,5));
		panel.setPreferredSize(new Dimension(190,80));
		
		uniformPrior = new JRadioButton(Strings.getString("uniformPrior"));
		uniformPrior.setPreferredSize(new Dimension(60,20));
		uniformPrior.setBorder(BorderFactory.createEmptyBorder(3,0,3,0));		
		uniformPrior.addActionListener(this);
		
		singlePrior = new JRadioButton(Strings.getString("singlePrior"));
		singlePrior.setPreferredSize(new Dimension(60,20));
		singlePrior.setBorder(BorderFactory.createEmptyBorder(3,0,3,0));
		singlePrior.addActionListener(this);
		radioButtonGroup = new ButtonGroup();
		radioButtonGroup.add(uniformPrior);
		radioButtonGroup.add(singlePrior);
		
		
		if (Parameters.prior_dirs.length == 0) {
			priorFilePanel = new FilePanel(".", true);
			uniformPrior.setSelected(true);
			priorFilePanel.setEnabled(false);
		}
		else {
			priorFilePanel = new FilePanel(Parameters.prior_dirs[0], true);			
			singlePrior.setSelected(true);
		}
		
		panel.add(uniformPrior);
		panel.add(singlePrior);
		panel.add(priorFilePanel);		
		this.add(panel, BorderLayout.CENTER);
	}

	
	/** Implements the actionPerformed function (from ActionListener) */
	public void actionPerformed(ActionEvent ev) {
		if (ev.getSource() == singlePrior) {
			priorFilePanel.setEnabled(true);
 			return;
 		}
		if (ev.getSource() == uniformPrior) {
			singlePrior.setEnabled(true);
			priorFilePanel.setEnabled(false);
 			return;
 		}
	}
}
