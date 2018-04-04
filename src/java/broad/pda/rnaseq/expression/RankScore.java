package broad.pda.rnaseq.expression;

import java.io.IOException;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;

public class RankScore {
	
	
	public RankScore(MatrixWithHeaders expression, String save)throws IOException{
		MatrixWithHeaders ranked=expression.rank();
		ranked.writeGCT(save);
	}

	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>1){
			MatrixWithHeaders expression=new MatrixWithHeaders(args[0]);
			String save=args[1];
			new RankScore(expression, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=gct file \n args[1]=save file";
	
}
