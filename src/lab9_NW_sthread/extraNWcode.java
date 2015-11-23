
/*
	
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
		
		System.out.println("---------------");
		//check traceback contents
		for( int x=0; x<tracebackMatrix.length; x++ )			//row index
		{
			String rows = "";
			for( int y=0; y<tracebackMatrix[0].length; y++ )		//col index
			{
				rows = rows + " " + tracebackMatrix[x][y];
			}
			System.out.println(rows);
		}
		
		
		
		
		
		
		
		
		
		
				//potential other 'best' alignments possible
				else if( maxScore == topScore && maxScore == diaScore )
				{
					tracebackMatrix[y-1][x-1] = "td";
				}
				else if( maxScore == topScore && maxScore == leftScore )
				{
					tracebackMatrix[y-1][x-1] = "tl";
				}
				else if( maxScore == diaScore && maxScore == leftScore )
				{
					tracebackMatrix[y-1][x-1] = "dl";
				}
				else if( maxScore == topScore && maxScore == diaScore && maxScore == leftScore )
				{
					tracebackMatrix[y-1][x-1] = "tdl";
				} 
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
*/
