package priority.gui;
import java.awt.*;
import javax.swing.*;
import javax.swing.text.StyledDocument;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.BadLocationException;

		
/**
 * The text panel.
 * @author raluca
 */
class TextPanel extends JPanel 
{
	static final long serialVersionUID = 1;

	private JScrollPane scrollPane;
	private JTextPane textPane;
	
    public StyledDocument doc;
    private Style textStyle;

	private static Color color1;
	private static Color color2;
	private static Color defaultColor = Color.BLACK;
	private static boolean odd = true;

	
	/** Constructor */
	public TextPanel() {
		super();
	}
	
	/** Creates and initilizes the components. */
	public void init() {
		this.setLayout( new BorderLayout());
		this.setBackground(Color.WHITE);
		this.setForeground(Color.WHITE);
		
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
		textPane.setMargin(new Insets(5,5,5,5));
		textPane.setBackground(Color.WHITE);
		textPane.setEditable(false);
		doc = (StyledDocument)textPane.getDocument(); /* get the text pane's document */
		textStyle = doc.addStyle("TextStyle", null);  /* create a style object */
		
		scrollPane = new JScrollPane(textPane, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
		                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		this.add(scrollPane, BorderLayout.CENTER);
		(textPane.getParent()).setBackground(Color.WHITE);
		
		float[] hsb = new float[3];
		Color.RGBtoHSB(180,0,0,hsb);
		color1 = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
		Color.RGBtoHSB(0,0,180,hsb);
		color2 = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
	}
	
	/** Returns the text in the text area. */
	public String getText() {
		return textPane.getText();
	}
	
	/** Appends a string to the text in the text area. */
	synchronized public void appendText(String str) 
	{
		try {
	        if (str.startsWith("-----")) /* should change in next release */
	        	odd = !odd;	    
	        
	        /* set the style color attribute */
	        if (odd)
	        	StyleConstants.setForeground(textStyle, color1);
	        else
	        	StyleConstants.setForeground(textStyle, color2);
	        
	        /* append the text to the document */
	        doc.insertString(doc.getLength(), str, textStyle);
	    } 
		catch (BadLocationException e) {
	    }
		/* make the last line visible */
	    textPane.setCaretPosition(textPane.getDocument().getLength()); 
	}

	/** Appends a string to the text in the text area, using the specified color. */
	synchronized public void appendText(String str, Color color) 
	{
		if (color == null)
			color = defaultColor;
		try {
			Color prev = StyleConstants.getForeground(textStyle);
        	StyleConstants.setForeground(textStyle, color);
	        
	        /* append the text to the document */
	        doc.insertString(doc.getLength(), str, textStyle);
	        
	        StyleConstants.setForeground(textStyle, prev);
	    } 
		catch (BadLocationException e) {
	    }
		/* make the last line visible */
	    textPane.setCaretPosition(textPane.getDocument().getLength()); 
	}
	
	/** Resets the text area. */
	synchronized public void resetText() {
		textPane.setText("");	
	}
}
