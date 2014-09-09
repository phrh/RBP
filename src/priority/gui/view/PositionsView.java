package priority.gui.view;
import java.awt.Color;
import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.BorderFactory;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.ButtonGroup;
import javax.swing.text.StyledDocument;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.BadLocationException;
import priority.Priority;
import priority.Parameters;
import priority.Strings;

/**
* A panel that shows the occurrences of the final motif.
* @author raluca
*/
public class PositionsView extends JFrame implements ActionListener {
	
	static final long serialVersionUID = 1;

	private JScrollPane scrollPane;
	private JTextPane textPane;
	
	private JRadioButton wrapButton, nowrapButton;
	private ButtonGroup radioButtonGroup;
	
    private StyledDocument doc;
    private Style textStyle;

	private static Color color1, color2, bkgrcolor1; //, bkgrcolor2;
	private int bestZ[];
	private String seq[];
	private String seq_name[];
	
	
	/** Constructor */
	public PositionsView(String title)
	{
		super(title);
	}
	
	
	/** Creates and initilizes the components. */
	public void init(int bestZ[], String seq[], String seq_name[]) 
	{
		this.seq = seq;
		this.seq_name = seq_name;
		this.bestZ = bestZ;
		
		
		this.getContentPane().setLayout( new BorderLayout());
		this.setBackground(Color.WHITE);
		this.setForeground(Color.WHITE);
		this.setIconImage(Priority.icon);
		this.setSize(540, 370);
		
		float[] hsb = new float[3];
		Color.RGBtoHSB(0,0,180,hsb);
		color1 = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
		Color.RGBtoHSB(180,0,0,hsb);
		color2 = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
		Color.RGBtoHSB(210,210,210,hsb);
		bkgrcolor1 = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
		//bkgrcolor2 = Color.getHSBColor(hsb[0], hsb[1], hsb[2]);
		
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
		textPane.setMargin(new Insets(10,7,10,7));
		textPane.setBackground(Color.WHITE);
		textPane.setEditable(false);
		
		doc = (StyledDocument)textPane.getDocument(); /* get the text pane's document */
		textStyle = doc.addStyle("TextStyle", null);  /* create a style object */
				
		scrollPane = new JScrollPane(textPane, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
		                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		getContentPane().add(scrollPane, BorderLayout.CENTER);
		(textPane.getParent()).setBackground(Color.WHITE);
		
		JPanel radioPanel = new JPanel(new GridLayout(1,2,10,10));
		radioPanel.setBorder(BorderFactory.createEmptyBorder(0,100,0,100));
		wrapButton = new JRadioButton(Strings.getString("wrapButton"));
		wrapButton.setPreferredSize(new Dimension(60,25));
		wrapButton.setBorder(BorderFactory.createEmptyBorder(5,30,5,0));		
		wrapButton.addActionListener(this);
		
		nowrapButton = new JRadioButton(Strings.getString("nowrapButton"));
		nowrapButton.setPreferredSize(new Dimension(60,25));
		nowrapButton.setSelected(true);
		nowrapButton.setBorder(BorderFactory.createEmptyBorder(5,0,5,30));
		nowrapButton.addActionListener(this);
		radioButtonGroup = new ButtonGroup();
		radioButtonGroup.add(wrapButton);
		radioButtonGroup.add(nowrapButton);
		
		radioPanel.add(wrapButton);
		radioPanel.add(nowrapButton);
		getContentPane().add(radioPanel, BorderLayout.NORTH);
		
		viewPositions(bestZ, seq, seq_name );
	}


	/** Appends a string to the text in the text area. */
	public void appendText(String str) 
	{
  		try 
  		{        
	        /* append the text to the document */
	        doc.insertString(doc.getLength(), str, textStyle);
	    } 
		catch (BadLocationException e) {}
	}

	
	/** Prints the occurences of the motif in the text area. */
	public void viewPositions(int bestZ[], String seq[], String seq_name[]) 
	{
		Color color;
		int position;
		for (int i=0; i<seq.length; i++) 
		{
			StyleConstants.setForeground(textStyle, Color.BLACK);
			appendText(seq_name[i] + " - ");
			
			position = bestZ[i];
		
			if (bestZ[i] < 0) 
			{
				appendText("no occurrence of the motif\n");
				appendText(seq[i] + "\n\n");
				continue;
			}
			else 
			{	
				if (bestZ[i] < seq[i].length()) 
				{/* direct strand */
				/* +1 because the index in a string starts from 0, but in the sequence it starts from 1 */
				appendText("the motif starts at position " + (bestZ[i]+1) + ", direct strand\n");
				color = color1;
				}
				else 
				{ /* complementary strand */
					appendText("the motif starts at position " + (bestZ[i]-seq[i].length()+1) + ", complementary strand\n");
					position = 2*seq[i].length() - bestZ[i] - Parameters.wsize;
					color = color2;
				}
			}
				
			appendText(seq[i].substring(0,position));
			StyleConstants.setBackground(textStyle, bkgrcolor1);
			StyleConstants.setForeground(textStyle, color);
			//System.out.println(seq[i].length() + " " + bestZ[i] + " " + position + " " + (position + Parameters.wsize));
			appendText( seq[i].substring(position, position + Parameters.wsize));
			
			StyleConstants.setBackground(textStyle, Color.WHITE);
			StyleConstants.setForeground(textStyle, Color.BLACK);
			appendText( seq[i].substring(position + Parameters.wsize, seq[i].length()) + "\n\n");
		}	
		
		/* make the first line visible */
	    textPane.setCaretPosition(0); 
	}

	
	/* To be modified in the next PRIORITY version */
	public void actionPerformed(ActionEvent ev) 
	{
		if (ev.getSource() == wrapButton) 
		{
			textPane = new JTextPane();
			textPane.setEditable(false);
			textPane.setFont(new Font("MonoSpaced", Font.PLAIN, 14));
			textPane.setMargin(new Insets(10,7,10,7));
			textPane.setBackground(Color.WHITE);
			textPane.setEditable(false);
			doc = (StyledDocument)textPane.getDocument(); /* get the text pane's document */
			textStyle = doc.addStyle("TextStyle", null);  /* create a style object */
					
			getContentPane().remove(scrollPane);
			scrollPane = new JScrollPane(textPane, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
			                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
			getContentPane().add(scrollPane, BorderLayout.CENTER);
			(textPane.getParent()).setBackground(Color.WHITE);
			viewPositions(bestZ, seq, seq_name );
			this.validate();
		}
		else 
		{ //nowrapButton
			textPane = new JTextPane()
			{
				private static final long serialVersionUID = 1;
				/* to disable the line wrapping in the JTextPane component */
				public void setSize(Dimension d)
				{
					if (d.width < getParent().getSize().width)
					{
						d.width = getParent().getSize().width;
					}
			
					super.setSize(d);
				}
				public boolean getScrollableTracksViewportWidth()
				{
					return false;
				}
			};
		
			textPane.setEditable(false);
			textPane.setFont(new Font("MonoSpaced", Font.PLAIN, 14));
			textPane.setMargin(new Insets(10,7,10,7));
			textPane.setBackground(Color.WHITE);
			textPane.setEditable(false);
			
			doc = (StyledDocument)textPane.getDocument(); /* get the text pane's document */
			textStyle = doc.addStyle("TextStyle", null);  /* create a style object */
			getContentPane().remove(scrollPane);
			scrollPane = new JScrollPane(textPane, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
			                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

			getContentPane().add(scrollPane, BorderLayout.CENTER);
			(textPane.getParent()).setBackground(Color.WHITE);
			viewPositions(bestZ, seq, seq_name );
			this.validate();
		}
	}
}