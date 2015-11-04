package lab7_AAquizGUI;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Random;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;

public class QuizShell extends JFrame
{
private static final long serialVersionUID = -3485284893489L;
	
	public JLabel question = new JLabel(" What is the code for:");
	public JLabel longName = new JLabel();
	public JTextField enterAnswer = new JTextField("");
	public JTextArea result = new JTextArea();
	private JButton cancelButton = new JButton("Cancel");
	private JButton submitButton = new JButton("Submit");
	public JLabel timer = new JLabel();
	
	private final IntroShell parent;
	private final QuizShell thisQuizShell;
	
	public static final Random RANDOM = new Random();
	private int rand = RANDOM.nextInt(20);
	public static final String[] SHORT_NAMES = {
		"A","R", "N", "D", "C", "Q", "E",
		"G",  "H", "I", "L", "K", "M", "F",
		"P", "S", "T", "W", "Y", "V" };
	public static final String[] FULL_NAMES = {
		"alanine","arginine", "asparagine",
		"aspartic acid", "cysteine",
		"glutamine",  "glutamic acid",
		"glycine" ,"histidine","isoleucine",
		"leucine",  "lysine", "methionine",
		"phenylalanine", "proline",
		"serine","threonine","tryptophan",
		"tyrosine", "valine" };
	HashMap<String, String> aa_Hash = new HashMap<>();
	{
		for( int x=0; x<SHORT_NAMES.length; x++)
		{
			aa_Hash.put(SHORT_NAMES[x], FULL_NAMES[x]);
		}
	}
	
	
	private void updateText()
	{
		result.setText(" The answer to the previous question will display here");
		result.setForeground(Color.BLACK);
		result.setEditable(false);
		longName.setText(FULL_NAMES[rand]+"?");
		
		validate();
	}
	
	private JPanel getCenterPanel()
	{
		JPanel tpanel = new JPanel();
		tpanel.setLayout(new GridLayout(1,3));
		tpanel.add(question);
		tpanel.add(longName);
		tpanel.add(enterAnswer);
		
		return tpanel;
	}
	private JPanel getBottomPanel()
	{
		JPanel bpanel = new JPanel();
		bpanel.setLayout(new GridLayout(1,3));
		bpanel.add(timer);
		bpanel.add(submitButton);
		bpanel.add(cancelButton);
		
		submitButton.addActionListener(new SubmitActionListener());
		submitButton.setMnemonic(KeyEvent.VK_ENTER);
		cancelButton.addActionListener(new CancelActionListener());
		
		return bpanel;
	}
	
	public QuizShell(IntroShell parent)
	{
		setTitle("Amino Acid Quiz: Goodluck!");
		this.parent = parent;
		this.thisQuizShell = this;
		setLocationRelativeTo(null);
		setSize(500,150);
		
		//content
		getContentPane().setLayout(new BorderLayout());
		getContentPane().add(result, BorderLayout.NORTH);
		getContentPane().add(getCenterPanel(), BorderLayout.CENTER);
		getContentPane().add(getBottomPanel(), BorderLayout.SOUTH);
		
		if( tr != null )
		{
			tr.setCancel();
		}
		tr = new TimerRunnable();
		new Thread(tr).start();
		
		
		updateText();
		setVisible(true);
		//force VM to kill frame in memory
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
//RUNNABLES************************************************************************************************
	protected TimerRunnable tr = null;
	public class TimerRunnable implements Runnable
	{
		private volatile boolean cancel = false;
		
		public void setCancel()
		{
			cancel = true;
		}
		
		public void run()
		{
			int time =30;
			
			try
			{
				time = Integer.parseInt(parent.enterTime.getText());
			}
			catch(Exception e)
			{
				e.printStackTrace();
			}
			
			while( time>=0 && ! cancel )
			{
				try
				{
					setTimeOnLabel(time);
					Thread.sleep(1000);
				}
				catch(Exception ex)
				{
					ex.printStackTrace();
				}
				time--;
			}
		}
	}
	
	private void setTimeOnLabel(final int time) throws Exception
	{
		SwingUtilities.invokeLater(new Runnable()
		{
			
			@Override
			public void run()
			{	
				if( time == 0 )
				{
					new ReportShell(thisQuizShell);
					dispose();
				}
				else
				{
					timer.setText(time + " seconds remaining");
				}
			}
		});
	}
	
//COUNTERS
	public int rightCount = 0;
	public int wrongCount = 0;
	public int total = rightCount+wrongCount;
	//hashmap for storing answers
	public LinkedHashMap<String, String> answer_Hash = new LinkedHashMap<>();
		
/*	public class QuizRunnable implements Runnable
	{
		public void run()
		{
			try
			{
				if( FULL_NAMES[rand].equals(aa_Hash.get(enterAnswer.getText().toUpperCase())))
				{
					rightCount++;
					total = rightCount+wrongCount;
					result.setText("  Correct! The code for "+FULL_NAMES[rand]+" is "+enterAnswer.getText().toUpperCase()+".");
					result.setForeground(Color.GREEN);
					answer_Hash.put("\nQuestion#"+total, "What is the code for: "+FULL_NAMES[rand]+"?\n	You answered:	"+enterAnswer.getText().toUpperCase()+"\n	Result:	Correct!");
					rand = RANDOM.nextInt(20);
					longName.setText(FULL_NAMES[rand]);
				}
				else
				{
					wrongCount++;
					total = rightCount+wrongCount;
					result.setText("  WRONG! The single letter code for "+FULL_NAMES[rand]+" is "+SHORT_NAMES[rand]+".");
					result.setForeground(Color.RED);
					answer_Hash.put("\nQuestion#"+total, "What is the code for: "+FULL_NAMES[rand]+"?\n	You answered:	"+enterAnswer.getText().toUpperCase()+"\n	Result:	WRONG! The code for "+FULL_NAMES[rand]+" is "+SHORT_NAMES[rand]);
					rand = RANDOM.nextInt(20);
					longName.setText(FULL_NAMES[rand]);
				}
				longName.setText(FULL_NAMES[rand]);
				enterAnswer.setText("");
				enterAnswer.requestFocus();
			}
			catch(Exception ex)
			{
				ex.printStackTrace();
			}
		}
	}
	*/
	
	//ACTION LISTENERS*****************************************************************************************


	public class SubmitActionListener implements ActionListener
	{
		public void actionPerformed(ActionEvent e)
		{
			//new Thread(new QuizRunnable()).start();
			
			if( FULL_NAMES[rand].equals(aa_Hash.get(enterAnswer.getText().toUpperCase())))
			{
				rightCount++;
				total = rightCount+wrongCount;
				result.setText("  Correct! The code for "+FULL_NAMES[rand]+" is "+enterAnswer.getText().toUpperCase()+".");
				result.setForeground(Color.GREEN);
				answer_Hash.put("\nQuestion#"+total, "What is the code for: "+FULL_NAMES[rand]+"?\n	You answered:	"+enterAnswer.getText().toUpperCase()+"\n	Result:	Correct!");
				rand = RANDOM.nextInt(20);
				longName.setText(FULL_NAMES[rand]);
			}
			else
			{
				wrongCount++;
				total = rightCount+wrongCount;
				result.setText("  WRONG! The single letter code for "+FULL_NAMES[rand]+" is "+SHORT_NAMES[rand]+".");
				result.setForeground(Color.RED);
				answer_Hash.put("\nQuestion#"+total, "What is the code for: "+FULL_NAMES[rand]+"?\n	You answered:	"+enterAnswer.getText().toUpperCase()+"\n	Result:	WRONG! The code for "+FULL_NAMES[rand]+" is "+SHORT_NAMES[rand]);
				rand = RANDOM.nextInt(20);
				longName.setText(FULL_NAMES[rand]);
			}
			longName.setText(FULL_NAMES[rand]);
			enterAnswer.setText("");
			enterAnswer.requestFocus();
		}
	}

	public class CancelActionListener implements ActionListener
	{
		public void actionPerformed(ActionEvent e)
		{
			if( tr != null )
			{
				tr.setCancel();
			}
			new ReportShell(thisQuizShell);
			dispose();
		}
	}
}
