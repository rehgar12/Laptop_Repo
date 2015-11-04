package lab7_AAquizGUI;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

public class ReportShell extends JFrame
{
	private static final long serialVersionUID = 6543345678998L;
	
	private final QuizShell parent;
	
	private JLabel header = new JLabel("How did you do?");
	private JTextArea report = new JTextArea();
	JScrollPane scroll = new JScrollPane (report, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
	private JButton redoButton = new JButton("Try Again");	//re-runs quizWindow with previously set time
	private static JButton closeButton = new JButton("Close");
	
	private void updateText()
	{
		report.setText("Summary Report: You answered "+parent.rightCount+" correct out of "+parent.total+" attempted.\n");
		//REPORT
		Object[] keyArray = parent.answer_Hash.keySet().toArray();
		for( Object key : keyArray )
		{
			report.append(key+"	"+parent.answer_Hash.get(key));
		}
		report.setEditable(false);
		validate();
	}
	
	private JPanel getBottomPanel()
	{
		JPanel bpanel = new JPanel();
		bpanel.setLayout(new GridLayout(0,2));
		bpanel.add(redoButton);	//tooltip, click or press Enter to re-start the quiz. 
		bpanel.add(closeButton);
		
		redoButton.addActionListener(new RedoActionListener());
		closeButton.addActionListener(new CloseActionListener());
		
		return bpanel;
	}
	
	public ReportShell(QuizShell parent)
	{
		this.parent = parent;
		setTitle("Amino Acid Quiz: Report");
		setLocationRelativeTo(null);
		setSize(500,450);
		
		//content
		getContentPane().setLayout(new BorderLayout());
		getContentPane().add(header, BorderLayout.NORTH);
		getContentPane().add(scroll, BorderLayout.CENTER);
		getContentPane().add(getBottomPanel(), BorderLayout.SOUTH);
		
		
		updateText();
		setVisible(true);
		//force VM to kill frame in memory
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
//RUNNABLES************************************************************************************************
/*	public class RedoRunnable implements Runnable
	{
		public void run()
		{
			try
			{
				new IntroShell();
				IntroShell.enterTime.requestFocus();
				dispose();
				//variable resets
				//Intro time
				IntroShell.enterTime.setText("");
				//Quiz counters and report hashmap
				parent.rightCount = 0;
				parent.wrongCount = 0;
				parent.total = 0;
				parent.answer_Hash.clear();
			}
			catch(Exception ex)
			{
				ex.printStackTrace();
			}
		}
	}*/
	
//ACTION LISTENERS*****************************************************************************************
	public class RedoActionListener implements ActionListener
	{
		public void actionPerformed(ActionEvent e)
		{
			new IntroShell();
			dispose();
		}
	}
	public class CloseActionListener implements ActionListener
	{
		public void actionPerformed(ActionEvent e)
		{
			System.exit(0);
		}
	}	
}

