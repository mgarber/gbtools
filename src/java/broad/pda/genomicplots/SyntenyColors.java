package broad.pda.genomicplots;

import java.awt.Color;
import java.lang.reflect.Field;

public class SyntenyColors {

	
	public static final Color CX = new Color(224,225,214);
	
	public static final Color CUn = Color.DARK_GRAY;
	
	public static final Color C1 = new Color(139, 155, 187);
	
	public static final Color CI = new Color(139, 155, 187);
	
	public static final Color C2 = new Color(206, 61, 50);
	
	public static final Color CII = new Color(206, 61, 50);
	
	public static final Color C2a =  C2.brighter();
	
	public static final Color C2b = C2a.brighter();
	
	public static final Color C3 = new Color(116, 155, 88);
	
	public static final Color CIII = new Color(116, 155, 88);
	
	public static final Color C4 = new Color(240, 230, 133);
	
	public static final Color CIV = new Color(240, 230, 133);
	
	public static final Color C5 = new Color(70, 105, 131);
	
	public static final Color C6 = new Color(186, 99, 56);
	
	public static final Color C7 = new Color(93, 177, 221);

	public static final Color C8 = new Color(128, 34, 104);
	
	public static final Color C9 = new Color(107, 215, 107);
	
	public static final Color C10 = new Color(213, 149, 167);
	
	public static final Color C11 = new Color(146, 72, 34);
	
	public static final Color C12 = new Color(131, 123, 141);
	
	public static final Color C13 = new Color(199, 81, 39);
	
	public static final Color C14 = new Color(213, 143, 92);
	
	public static final Color C15 = new Color(122, 101, 165);
	
	public static final Color C16 = new Color(228, 175, 105);
	
	public static final Color C17 = new Color(59, 27, 83);
	
	public static final Color C18 = new Color(205, 222, 183);
	
	public static final Color C19 = new Color(97, 42, 121);
	
	public static final Color C20 = new Color(174, 31, 99);
	
	public static final Color C21 = new Color(231, 199, 111);
	
	public static final Color C22 = new Color(90, 101, 94);
	
	public static final Color C23 = new Color(204, 153, 0);
	
	public static final Color C24 = new Color(153, 204, 0);
	
	public static final Color C25 = new Color(51, 204, 0);
	
	public static final Color C26 = new Color(0, 204, 51);
	
	public static final Color C27 = new Color(0, 204, 153);
	
	public static final Color C28 = new Color(0, 153, 204);
	
	public static final Color C29 = new Color(10, 71, 255);
	
	public static final Color C30 = new Color(71, 117, 255);
	
	public static final Color C31 = new Color(255, 194, 10);
	
	public static final Color C32 = new Color(255, 209, 71);
	
	public static final Color C33 = new Color(153, 0, 51);
	
	public static final Color C34 = new Color(153, 26, 0);
	
	public static final Color C35 = new Color(153, 102, 0);
	
	public static final Color C36 = new Color(128, 153, 0);
	
	public static final Color C37 = new Color(51, 153, 0);
	
	public static final Color C38 = new Color(0, 153, 26);
	
	public static final Color C39 = new Color(0, 153, 102);
	
	public static final Color C40 = new Color(0, 128, 153);
	
	public static final Color C41 = new Color(0, 51, 153);
	
	public static final Color C42 = new Color(26, 0, 153);
	
	public static final Color C43 = new Color(102, 0, 153);
	
	public static final Color C44 = new Color(153, 0, 128);
	
	public static final Color C45 = new Color(214, 0, 71);
	
	public static final Color C46 = new Color(255, 20, 99);
	
	public static final Color C47 = new Color(0, 214, 143);
	
	public static final Color C48 = new Color(20, 255, 177);
	
	public static final Color darkGreen = new Color(0, 153, 0);
	
	public static final Color violet = new Color(204, 0, 255);
	
	public static final Color lightGreen = new Color(102, 255, 0);
	
	public static final Color brown = new Color(102, 51, 0);
	
	public static final Color paleDullOrange = new Color(255, 204, 153);
	
	public static final Color class1 = Color.red;
	
	public static final Color class2 = darkGreen;
	
	public static final Color class3 = Color.cyan;
	
	//public static final Color class4 = violet;
	
	//public static final Color class5 = lightGreen;
	
	//public static final Color class6 = Color.blue;
	
	//public static final Color class7 = Color.orange;
	
	//public static final Color class8 = brown;
	
	//public static final Color class9 = paleDullOrange;
	
	//public static final Color class10 = Color.pink;
	
	//public static final Color class11 = Color.magenta;
	
	public static final Color dogBreak = new Color(249, 210, 249); //darker //Color(255,230,255); //lighter 
	
	public static final Color humanBreak = new Color(230,230,255); //ligher (193, 191, 221); //darker
	
	public static final Color mouseBreak = new Color(255,255,230); //lighter (198, 247, 246); //darker
	
	public static final Color rodentBreak = new Color(230,255,255); //lighter // (187, 249, 187); darker
	
	public static final Color chimpBreak = new Color(255, 230, 230);
	
	public static final Color ratBreak = new Color(230,255,230); //darker (230, 255, 230); //lighter
	
	public static Color getChromosomeColor(String chromosomeNumber) {
		Class c;
		Color result = Color.BLACK;
		
		try {
			c = Class.forName("edu.mit.broad.prodinfo.genomicplot.SyntenyColors");
			Field f = c.getField("C" + chromosomeNumber);
			result = (Color) f.get(new SyntenyColors());
		} catch (Exception e) {
			//System.err.println("Chromosome number " + chromosomeNumber + " has no color defined returning default");
			result = Color.LIGHT_GRAY;
		}

		
		return result;
	}
	
	public static Color getLinageColor(String linage)
	throws ClassNotFoundException, SecurityException,
	NoSuchFieldException, IllegalArgumentException,
	IllegalAccessException {
		Class c;
		Color result = Color.BLACK;
		
		c = Class.forName("broad.pda.genomicplot.SyntenyColors");
		Field f = c.getField(linage + "Break");
		result = (Color) f.get(new SyntenyColors());
		
		return result;
	}
	
	public static Color getClassColor(String dupClass)
	throws ClassNotFoundException, SecurityException,
	NoSuchFieldException, IllegalArgumentException,
	IllegalAccessException {
		Class c;
		Color result = Color.BLACK;
		c = Class.forName("broad.pda.genomicplots.SyntenyColors");
		Field f = c.getField(dupClass);
		result = (Color) f.get(new SyntenyColors());
		
		return result;
	}
	
}