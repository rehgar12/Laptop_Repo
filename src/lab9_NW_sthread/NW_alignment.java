package lab9_NW_sthread;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class NW_alignment
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
		String[] seq1Array = seq1.split("");	//[0] is ""
		String[] seq2Array = seq2.split("");	//[0] is ""
		int seqMatrix[][] = new int[seq1Array.length][seq2Array.length];
		String tracebackMatrix[][] = new String[seq1Array.length-1][seq2Array.length-1];
		//initialize and fill in using NW and relevant scoring matrix
		int colInit = -1;
		int rowInit = 0;
		
		long startNW = System.currentTimeMillis();
		for( int y=0; y<seqMatrix.length; y++ )			//row index (y-axis)
		{
			for( int x=0; x<seqMatrix[0].length; x++ )	//col index (x-axis)
			{
				if( y == 0 )		//first row initializer
				{
					seqMatrix[y][x] = rowInit;
//					rowInit--;
					rowInit = rowInit - 8;
				}
				else
				{
					if( x == 0)		//first col initializer
					{
						seqMatrix[y][x] = colInit;
//						colInit--;
						colInit = colInit - 8;
					}
					else		//fill in most positive score
					{
						String aaPair = seq1Array[y] + seq2Array[x];	//get value for pair[y][x] from scoring matrix
						int alignScore = Integer.parseInt(scoreMat.get(aaPair));
						int topScore = 0;			//top		y-1
						topScore = alignScore - 8 + seqMatrix[y-1][x];							
						int diaScore = 0;			//diagonal	y-1, x-1
						diaScore = alignScore + seqMatrix[y-1][x-1];
						int leftScore = 0;			//left		x-1
						leftScore = alignScore - 8 + seqMatrix[y][x-1];
						int maxScore = Math.max(topScore, Math.max(diaScore, leftScore));	//set value of [y][x] to maxScore
						seqMatrix[y][x] = maxScore;
						if( maxScore == topScore )			//populate traceback matrix
						{
							tracebackMatrix[y-1][x-1] = "t";	//space in seq2
						}
						else if( maxScore == diaScore )
						{
							tracebackMatrix[y-1][x-1] = "d";	//no space
						}
						else if( maxScore == leftScore )
						{
							tracebackMatrix[y-1][x-1] = "l";	//space in seq1
						}
					}
				}
			}
		}
		long endNW = System.currentTimeMillis();
		System.out.println("NW part took " + (endNW - startNW) + " milliseconds.");
		
		long startTB = System.currentTimeMillis();
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
				if( tracebackMatrix[x][y].equals("t"))
				{
					seq1Align = seq1Align + seq1Array[x+1];
					seq2Align = seq2Align + "-";
				}
				else if( tracebackMatrix[x][y].equals("d"))
				{
					seq1Align = seq1Align + seq1Array[x+1];
					seq2Align = seq2Align + seq2Array[y+1];
					colCount++;
				}
				else if( tracebackMatrix[x][y].equals("l"))
				{
					seq1Align = seq1Align + "-";
					seq2Align = seq2Align + seq2Array[y+1];
					colCount++;
					x++;		//reset x to stay in the same row
				}
			}
		}
		long endTB = System.currentTimeMillis();
		System.out.println("TB took " + (endTB - startTB) + " milliseconds.");
		
		
		String seq1Alignment = new StringBuffer(seq1Align).reverse().toString();
		String seq2Alignment = new StringBuffer(seq2Align).reverse().toString();
		System.out.println(seq1Alignment);
		System.out.println(seq2Alignment);
	}
	
	public static void main(String[] args) throws Exception
	{
		long startTime = System.currentTimeMillis();
		//get scoring matrix
		HashMap<String, String> scoreMat = getScoringMatrix("/home/playerra/Documents/UNCC_fall2015/binf6380/lab9/blosum50.mtx");
		//get seqs to align
		HashMap<String, String> seqHash = getSeqs("/home/playerra/Documents/UNCC_fall2015/binf6380/lab9/twoSeqs.fasta");
		//get seqs from fasta or copy/paste into seq1 and seq2 values
		String seq1 = seqHash.get("seq1");
		String seq2 = seqHash.get("seq2");	
		System.out.println("Aligning using single-threaded Needleman-Wunsch algorithm...\n");
		//the alignment
//		long startTime = System.currentTimeMillis();
		NW_algorithm(scoreMat, seq1, seq2);
		long endTime = System.currentTimeMillis();
		System.out.println("\nThis alignment took " + (endTime - startTime) + " milliseconds.");
	}
}