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
					//System.out.println(headerItem + "|" + rowItem);
					scoreList.add(rowItem);
					checkIndexCount++;
				}
			}	
			
			if( checkIndexCount == 20 )
			{
				int itemCount = 0;
				for( String scoreItem : scoreList )
				{
					//System.out.println(headersList.get(headerCount) + headersList.get(itemCount)+ scoreItem);
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
		//get fasta sequences from infile
		BufferedReader reader = new BufferedReader(new FileReader(inFile));
		String seq1 = "";
		String seq2 = "";
		//skipping first line, the following for loops restrict the alignment to 2 seqs
		reader.readLine();
		for (String line = reader.readLine(); !line.startsWith(">"); line = reader.readLine()) 
		{
			seq1 = seq1 + line;
		}
		for (String line = reader.readLine(); line != null; line = reader.readLine()) 
		{
			seq2 = seq2 + line;
		}
		reader.close();
		seqHash.put("seq1", seq1);
		seqHash.put("seq2", seq2);
		return seqHash;
	}
	
	
	public static void NW_algorithm(HashMap<String, String> scoreMat, HashMap<String, String> seqHash)
	{
		String seq1 = "ACDEFG";
		String seq2 = "HIKL";
		String[] seq1Array = seq1.split("");	//[0] is ""
		String[] seq2Array = seq2.split("");	//[0] is ""
		
		int seqMatrix[][] = new int[seq1Array.length][seq2Array.length];
		String tracebackMatrix[][] = new String[seq1Array.length][seq2Array.length];
		//print num rows and num cols
		System.out.println(seqMatrix.length + " " + seqMatrix[0].length);
		
	//	System.out.println(seqMatrix[0][1]);
		
		//initialize
		int colInit = -1;
		int rowInit = 0;
		for( int y=0; y<seqMatrix.length; y++ )			//row index (y-axis)
		{
			for( int x=0; x<seqMatrix[0].length; x++ )	//col index (x-axis)
			{
				if( y == 0 )		//first row initializer
				{
					seqMatrix[y][x] = rowInit;
					rowInit--;
				}
				else
				{
					if( x == 0)		//first col initializer
					{
						seqMatrix[y][x] = colInit;
						colInit--;
					}
					else		//fill in most positive score, track traceback
					{
						//get value for pair[y][x] from scoring matrix
						String aaPair = seq1Array[y] + seq2Array[x];
						int alignScore = Integer.parseInt(scoreMat.get(aaPair));
						
						//top		y-1
						int topScore = 0;
						topScore = alignScore + seqMatrix[y-1][x];
						
						//diagonal	y-1, x-1
						int diaScore = 0;
						diaScore = alignScore + seqMatrix[y-1][x-1];
						
						//left		x-1
						int leftScore = 0;
						leftScore = alignScore + seqMatrix[y][x-1];
						
						//set value of [y][x] to maxScore
						int maxScore = Math.max(topScore, Math.max(diaScore, leftScore));
						seqMatrix[y][x] = maxScore;
					}
				}
			}
		}
		
//		System.out.println(seqMatrix[6][4]);
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		//check matrix contents
		for( int x=0; x<seqMatrix.length; x++ )			//row index
		{
			String rows = "";
			for( int y=0; y<seqMatrix[0].length; y++ )		//col index
			{
				rows = rows + " " + seqMatrix[x][y];
			}
			System.out.println(rows);
		}
		
		System.out.println(seq1Array[6]);
		
	}
	
	
	
	
	
	public static void main(String[] args) throws Exception
	{
		//get scoring matrix
		HashMap<String, String> scoreMat = getScoringMatrix("/home/playerra/Documents/UNCC_fall2015/binf6380/lab9/blosum50.mtx");
	/*	for (Map.Entry<String,String> entry : scoreMat.entrySet()) {
			  String key = entry.getKey();
			  String value = entry.getValue();
			  System.out.println(key + "|" + value);
			}
	*/		
		//get seqs to align
		HashMap<String, String> seqHash = getSeqs("/home/playerra/Documents/UNCC_fall2015/binf6380/lab9/twoSeqs.fasta");

		NW_algorithm(scoreMat, seqHash);
	}
}
