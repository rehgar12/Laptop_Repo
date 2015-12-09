package lab9_NW_mthreaded_EC;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.CountDownLatch;

public class NW_multiThread
{
	public static HashMap<String, String> getScoringMatrix(String matrixFilePath) throws Exception
	{
		HashMap<String, String> pairScore = new HashMap<>();
		BufferedReader reader = new BufferedReader(new FileReader(matrixFilePath));
		String line = "";
		line = reader.readLine();
		String[] headersArray;
		headersArray = line.split("\\W+");		//[0] is empty
		List<String> headersList = new ArrayList<>();
		for( String x : headersArray )
		{
			headersList.add(x);
		}
		headersList.remove(0);
		int headerCount = 0;
		while( (line = reader.readLine()) != null )
		{
			String[] scoreArray;
			scoreArray = line.split(" ");
			List<String> scoreListALL = new ArrayList<>();
			List<String> scoreList = new ArrayList<>();
			for( String x : scoreArray )
			{
				scoreListALL.add(x);
			}
			scoreListALL.remove(0);		//leave only scores in the list
			int checkIndexCount = 0;
			for( String rowItem : scoreListALL)
			{
				if( !rowItem.equals("") )
				{
					scoreList.add(rowItem);
					checkIndexCount++;
				}
			}	
			if( checkIndexCount == 20 )
			{
				int itemCount = 0;
				for( String scoreItem : scoreList )
				{
					pairScore.put(headersList.get(headerCount) + headersList.get(itemCount), scoreItem);
					itemCount++;
				}
			}
			headerCount++;
		}		
		reader.close();
		return pairScore;
	}
	
	public static HashMap<String, String> getSeqs(String inFile) throws Exception
	{
		HashMap<String, String> seqHash = new HashMap<>();
		BufferedReader reader = new BufferedReader(new FileReader(inFile)); 		//get fasta sequences from infile
		String seq1 = "";
		String seq2 = "";
		reader.readLine();		//skipping first line, the following for loops restrict the alignment to 2 seqs
		for (String line = reader.readLine(); !line.startsWith(">"); line = reader.readLine()) 
		{
			seq1 = seq1 + line;
		}
		for (String line = reader.readLine(); line != null; line = reader.readLine()) 
		{
			seq2 = seq2 + line;
		}
		reader.close();
		if( seq1.length() >= seq2.length() )	//seq1 must be the longer of the 2 seqs
		{
			seqHash.put("seq1", seq1);
			seqHash.put("seq2", seq2);
		}
		else
		{
			seqHash.put("seq1", seq2);
			seqHash.put("seq2", seq1);
		}
		return seqHash;
	}
	
	public static void NW_algorithm(HashMap<String, String> scoreMat, String seq1, String seq2)
	{
		final String[] seq1Array = seq1.split("");	//[0] is ""
		final String[] seq2Array = seq2.split("");	//[0] is ""
		final int seqMatrix[][] = new int[seq1Array.length][seq2Array.length];
		final char tracebackMatrix[][] = new char[seq1Array.length-1][seq2Array.length-1];
		final HashMap<String, String> finalScoreMat = scoreMat;
		
		//setup CountDownLatch, must wait on threads to finish populating the traceback matrix before visualizing the alignment
		int numThreads = seqMatrix.length-1;
		final CountDownLatch cdl = new CountDownLatch(numThreads);
				
		//initialize and fill in using NW and relevant scoring matrix
		int colInit = -1;
		int rowInit = 0;
		for( int y=0; y<seqMatrix.length; y++ )			//row index (y-axis)
		{
			final int yThread = y;						//to be passed into thread
			for( int x=0; x<seqMatrix[0].length; x++ )	//col index (x-axis)
			{
				final int xThread = x;					//to be passed into thread
				if( y == 0 )		//first row initializer
				{
					seqMatrix[yThread][xThread] = rowInit;
//					rowInit--;
					rowInit = rowInit - 8;
//					System.out.println(seqMatrix[y][x]);
				}
				else
				{
					if( x == 0)		//first col initializer
					{
						seqMatrix[yThread][xThread] = colInit;
//						colInit--;
						colInit = colInit - 8;
//						System.out.println(seqMatrix[y][x]);
					}
					else	//spawn a thread to calc scores for each row (after first row, y=0, and first col, x=0, have been initialized)	
					//fill in most positive score
					{
						Thread rowThread = new Thread(new Runnable() 
						{
							public void run() 
							{
								try
								{
									
									String aaPair = seq1Array[yThread] + seq2Array[xThread];	//get value for pair[y][x] from scoring matrix
									int alignScore = Integer.parseInt(finalScoreMat.get(aaPair));
									int topScore = 0;			//top		y-1	(gap penalty of -8)
									topScore = alignScore - 8 + seqMatrix[yThread-1][xThread];							
									int diaScore = 0;			//diagonal	y-1, x-1
									diaScore = alignScore + seqMatrix[yThread-1][xThread-1];
									int leftScore = 0;			//left		x-1 (gap penalty of -8)
									leftScore = alignScore - 8 + seqMatrix[yThread][xThread-1];
									int maxScore = Math.max(topScore, Math.max(diaScore, leftScore));	//set value of [y][x] to maxScore
									seqMatrix[yThread][xThread] = maxScore;
									if( maxScore == topScore )			//populate traceback matrix
									{
										tracebackMatrix[yThread-1][xThread-1] = 't';	//space in seq2
									}
									else if( maxScore == diaScore )
									{
										tracebackMatrix[yThread-1][xThread-1] = 'd';	//no space
									}
									else if( maxScore == leftScore )
									{
										tracebackMatrix[yThread-1][xThread-1] = 'l';	//space in seq1
									}
									cdl.countDown();
								}
								catch(Exception e)
								{
									e.printStackTrace();
								}
							}
						} );
						rowThread.start();
					}
				}
			}
		}
		
		try
		{
			cdl.await();
		}
		catch(InterruptedException ie)
		{
			System.out.println("wubbawubba");
		}
		
		//visualize alignment of seqs
		String seq1Align = "";
		String seq2Align = "";
		int colCount = 0;
		for( int x=tracebackMatrix.length-1; x>=0; x-- )			//row index
		{
			int y=tracebackMatrix[0].length-1;
			y = y - colCount;										//col index
			if( y<0 )
			{
				y = 0;	
				seq1Align = seq1Align + seq1Array[x+1];
				seq2Align = seq2Align + "-";
				colCount++;
			}
			else
				{
				if( tracebackMatrix[x][y] == 't' )
				{
					seq1Align = seq1Align + seq1Array[x+1];
					seq2Align = seq2Align + "-";
				}
				else if( tracebackMatrix[x][y] == 'd')
				{
					seq1Align = seq1Align + seq1Array[x+1];
					seq2Align = seq2Align + seq2Array[y+1];
					colCount++;
				}
				else if( tracebackMatrix[x][y] == 'l')
				{
					seq1Align = seq1Align + "-";
					seq2Align = seq2Align + seq2Array[y+1];
					colCount++;
					x++;		//reset x to stay in the same row
				}
			}
		}
		String seq1Alignment = new StringBuffer(seq1Align).reverse().toString();
		String seq2Alignment = new StringBuffer(seq2Align).reverse().toString();
		System.out.println(seq1Alignment);
		System.out.println(seq2Alignment);
	}
	
	public static void main(String[] args) throws Exception
	{
		long startTime1 = System.currentTimeMillis();
		//get scoring matrix
		final HashMap<String, String> scoreMat = getScoringMatrix("/home/playerra/Documents/UNCC_fall2015/binf6380/lab9/blosum50.mtx");
		long endTime1 = System.currentTimeMillis();
		System.out.println("Getting Scoring Matrix took " + (endTime1 - startTime1) + " milliseconds.");
	
		
		long startTime2 = System.currentTimeMillis();
		//get seqs to align
		HashMap<String, String> seqHash = getSeqs("/home/playerra/Documents/UNCC_fall2015/binf6380/lab9/twoSeqs.fasta");
		//get seqs from fasta or copy/paste into seq1 and seq2 values
		String seq1 = seqHash.get("seq1");
		String seq2 = seqHash.get("seq2");	
		long endTime2 = System.currentTimeMillis();
		System.out.println("Getting sequences to align took " + (endTime2 - startTime2) + " milliseconds.");
		
		
//test seqs
//		String seq1 = "EEFFYYYFF";
//		String seq2 = "GYYYD";	
		
		System.out.println("\nAligning using single-threaded Needleman-Wunsch algorithm...\n");
		//the alignment
		long startTime3 = System.currentTimeMillis();
		NW_algorithm(scoreMat, seq1, seq2);
		long endTime3 = System.currentTimeMillis();
		System.out.println("\nThis alignment took " + (endTime3 - startTime3) + " milliseconds.");
	}
}