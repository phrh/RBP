package priority.gui.view;
import java.awt.Color;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Insets;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;
import javax.swing.JFrame;
import javax.swing.text.StyledDocument;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.BadLocationException;

import priority.Priority;
import priority.gibbs.GibbsStatic;

/**
* A panel that shows the final motif.
* @author raluca
*/
public class MotifView extends JFrame {
	static final long serialVersionUID = 1;

	private JScrollPane scrollPane;
	private JTextPane textPane;
	
    private StyledDocument doc;
    private Style textStyle;

	private static Color color;
	private static final Color[] DNAcolors = {Color.RED, Color.BLUE, Color.ORANGE, Color.GREEN};
	private static final double ln2 = Math.log(2);
	
	/** Constructor */
	public MotifView(String title) {
		super(title);
	}
	
	/** Creates and initilizes the components. */
	public void init(double phi[][]) {
		this.getContentPane().setLayout( new BorderLayout());
		this.setBackground(Color.WHITE);
		this.setForeground(Color.WHITE);
		this.setIconImage(Priority.icon);
		this.setSize(670, 300);
		
		float[] hsb = new float[3];
		Color.RGBtoHSB(0,0,180,hsb);
		color = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
		color = Color.BLACK;
		
		textPane = new JTextPane()
		{
			private static final long serialVersionUID = 1;
			/* to disable the line wrapping in the JTextPane component */
			public void setSize(Dimension d)
			{
				if (d.width < getParent().getSize().width)
					d.width = getParent().getSize().width;
				super.setSize(d);
			}
			public boolean getScrollableTracksViewportWidth()
			{
				return false;
			}
		};
		textPane.setEditable(false);
		textPane.setFont(new Font("MonoSpaced", Font.PLAIN, 14));
		textPane.setMargin(new Insets(10,5,10,5));
		textPane.setBackground(Color.WHITE);
		textPane.setEditable(false);
		doc = (StyledDocument)textPane.getDocument(); /* get the text pane's document */
		textStyle = doc.addStyle("TextStyle", null);  /* create a style object */
		StyleConstants.setForeground(textStyle, color);
		
		scrollPane = new JScrollPane(textPane, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
		                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		getContentPane().add(scrollPane, BorderLayout.CENTER);
		(textPane.getParent()).setBackground(Color.WHITE);
		
		viewPSSM(phi);
	}


	/** Appends a string to the text in the text area. */
	public void appendText(String str) 
	{
  		try {        
	        /* append the text to the document */
	        doc.insertString(doc.getLength(), str, textStyle);
	    } 
		catch (BadLocationException e) {}
		/* make the last line visible */
	    textPane.setCaretPosition(textPane.getDocument().getLength()); 
	}
	

	/** Appends a string to the text in the text area, using a specific color. */
	public void appendText(String str, Color newcolor) 
	{
  		try {
  			StyleConstants.setForeground(textStyle, newcolor);
	        /* append the text to the document */
	        doc.insertString(doc.getLength(), str, textStyle);
	        StyleConstants.setForeground(textStyle, color);
	    } 
		catch (BadLocationException e) {}
		/* make the last line visible */
	    textPane.setCaretPosition(textPane.getDocument().getLength()); 
	}

	
	/** Writes the PSSM in the text area. */
	public void viewPSSM(double[][] phi) 
	{
		int dec = 4;
		appendText(GibbsStatic.repeatChar(' ', 22));
		for (int j=0; j<phi[0].length; j++)
			appendText(GibbsStatic.formatInt(j+1, dec+2) + " ");
		appendText("\n");
		for (int i=0; i<phi.length; i++)
		{
			if (i != 1)
				appendText(GibbsStatic.repeatChar(' ', 18) + GibbsStatic.DNAchars[i] + "  ");
			else 
				appendText(GibbsStatic.repeatChar(' ', 11) + "PSSM:  " + GibbsStatic.DNAchars[i] + "  ");
			for (int j=0; j<phi[i].length; j++)
				appendText(GibbsStatic.formatDouble01(phi[i][j], dec) + " ");
			appendText("\n");
		}	
			
		appendText("\nMultilevel consensus: ");
		
		
		StyleConstants.setFontSize(textStyle, 18);
		/* now write the letters with the highest frequencies */
		String letter; 
		int max;
		double phi_temp[][] = new double[phi.length][phi[0].length];
		for (int i=0; i<phi.length; i++)
			for (int j=0; j<phi[0].length; j++)
				phi_temp[i][j] = phi[i][j];
		
		String line = "";
		for (int k=0; k<4; k++) {
			line = "";
			for (int j=0; j<phi_temp[0].length; j++) {
				max = 0;
				for (int i=1; i<4; i++)
					if (phi_temp[i][j] > phi_temp[max][j]) max = i;
				letter = "" + GibbsStatic.DNAchars[max];
				if (phi_temp[max][j] < 0.2)
					if (line.trim().length() == 0)
						line = line + "     ";
					else
						appendText("     ");
				else {
					if (line.trim().length() == 0) {
						if (k>0) {
							StyleConstants.setFontSize(textStyle, 14);
							appendText("                      ");
							StyleConstants.setFontSize(textStyle, 18);
							appendText(line);
						}
						line = ".";
					}
					if (phi_temp[max][j] < 0.5)
						appendText("  " + letter.toLowerCase() + "  ", DNAcolors[max]);
					else
						appendText("  " + letter.toUpperCase() + "  ", DNAcolors[max]);
				}
				phi_temp[max][j] = -1;
				//StyleConstants.setFontSize(textStyle, 14);
			}
			if (line.trim().length() != 0) appendText("\n");				
		}
		StyleConstants.setFontSize(textStyle, 14);
		
		double[] entr = computeInformation(phi);
		appendText("\nInformation (bits):  ");
		for (int j=0; j<entr.length; j++)
			appendText(GibbsStatic.formatDouble01(2-entr[j], dec) + " ");
	}
	
	/** Computes the information content. */
	public double[] computeInformation(double[][] phi) {
		int size = phi[0].length;
		double entr[] = new double[size];
		
		for (int i=0; i<size; i++) {
			entr[i] = 0;
			for (int j=0; j<phi.length; j++)
				if (phi[j][i] > 0)
					entr[i] -= phi[j][i] * Math.log(phi[j][i])/ln2;
		}
		return entr;
	}
}
