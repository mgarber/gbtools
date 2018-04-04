package broad.pda.graph.dotUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.Map;
import java.util.TreeMap;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT.BubbleEdge;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT.EdgeSourceType;

public class ExtractRegionsFromDOT {

	private static final String quote="\"";
	double alpha=.05;
	
	public ExtractRegionsFromDOT(File DOTFile, Alignments region, String save)throws IOException, ParseException{
		ChromosomeWithBubblesJGraphT graph=buildFromDOT(DOTFile);
		String regionName = region.getName();
		regionName = regionName.replace(":", "_");
		regionName = regionName.replace("-", "_");//This substitutions are necessary for Graphbiz to correctly process the file
		graph.writeGraph(save, region, true, regionName);
		
		
		//Map<Path, double[]> pairedEndPaths= ContinuousDataAlignmentModel.scorePaths(graph.getAllPaths(), graph.getLambda(), graph.getNumberOfMarkers(), graph.getNumberOfReads(), graph.getRPKMConstant(), graph.getLocalRate(), alpha);
		//Collection<RefSeqGene> genes=ContinuousDataAlignmentModel.filterPaths(pairedEndPaths, .05);
		//BEDFileParser.writeFullBED(save+".bed", graph.getGenePaths(1));
	}
	
	public static ChromosomeWithBubblesJGraphT buildFromDOT(File DOTFile) throws IOException, ParseException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(DOTFile)));
		ChromosomeWithBubblesJGraphT graph = buildFromDOT(reader);
		reader.close();
		return graph;
	}
	
	public static ChromosomeWithBubblesJGraphT buildFromDOT(BufferedReader reader)  throws IOException, ParseException{
		return buildFromDOT(reader, 0);
	}
	
	public static ChromosomeWithBubblesJGraphT buildFromDOT(BufferedReader reader, int minSpliceNum) throws IOException, ParseException{
		String name="";
		double lambda=0;
		double numMarkers=0;
		double numReads=0;
		Map<LightweightGenomicAnnotation, Double> nodes=new TreeMap<LightweightGenomicAnnotation, Double>();
		Map<Alignments, Double> localRateMap=new TreeMap<Alignments, Double>();
		Map<BubbleEdge, Double> edges=new TreeMap<BubbleEdge, Double>();
		
    	String nextLine;
    	while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
    		//if line contains ->, then check if its visible, if not just grab nodes, else grab nodes and edges
    		String[] tokens=nextLine.split("->");
    		if(tokens.length>1){
    			
    			boolean isHidden=false;
    			String style=getStyle(nextLine);
    			if(style.equalsIgnoreCase("invis")){isHidden=true;}
    			    			
    			//grab nodes
    			for(int i=0; i<tokens.length; i++){
    				Alignments nodeAlign=getNode(tokens[i]);
    				//nodes.add(nodeAlign);
    				if(!isHidden && i<tokens.length-1){
        				Alignments left=nodeAlign;
        				Alignments right=getNode(tokens[i+1]);
        				Alignments connection=getConnection(left, right);
        				double count=getCount(nextLine);
        				EdgeSourceType type=getType(nextLine);
        				//System.err.println("Edge between " + left.toUCSC() + " and " + right.toUCSC() + " type " + type + " count " + count + " is type equals to SPLICED? " + EdgeSourceType.PAIRED.equals(type) + " really? " + (type.ordinal() == EdgeSourceType.SPLICED.ordinal() )+ " minSpliceNum? " + minSpliceNum + " passes? " + (type.ordinal() != EdgeSourceType.SPLICED.ordinal() ||(type.ordinal() == EdgeSourceType.SPLICED.ordinal() && count > minSpliceNum )));
        				if(type.ordinal() != EdgeSourceType.SPLICED.ordinal() || (type.ordinal() == EdgeSourceType.SPLICED.ordinal() && count > minSpliceNum ) ) {
        					BubbleEdge edge=new BubbleEdge(connection, left, right, count, type);
        					if(EdgeSourceType.PAIRED.ordinal() == type.ordinal()) {
        						edge.setPairedEndCounts(count);
        					}
        					edges.put(edge, count);
        				} //else {
        					//System.err.println("Edge " + connection.toUCSC() + " type " + type + " did not pass filter spliceCount: " + count);
        				//}
        				//BubbleEdge edge=new BubbleEdge(connection, left, right, count, type);
        				//edges.put(edge, count);
    				}
    			}
    			
    		}
    		else if(nextLine.startsWith("digraph")){
    			NumberFormat nf = NumberFormat.getInstance();
    			name=nextLine.split(" ")[1].split("_")[0];
    			try{
    			lambda=new Double(nextLine.split(" ")[1].split("_")[1]);
    			numMarkers=nf.parse(nextLine.split(" ")[1].split("_")[2]).doubleValue();
    			numReads=new Double(nextLine.split(" ")[1].split("_")[3]);
    			}catch(ArrayIndexOutOfBoundsException ex){}
    		}
    		else if(nextLine.startsWith(quote)){
    			Alignments node=getNode(nextLine.split(" ")[0]);
    			double count=getNodeCount(nextLine); //TODO Merge the read and write functions
    			double localRate=getLocalRate(nextLine);
    			nodes.put(node, count);
    			localRateMap.put(node, localRate);
    			//System.out.println(node.toUCSC()+"\t"+count);
    		}
    	}
    	
    	ChromosomeWithBubblesJGraphT graph=new ChromosomeWithBubblesJGraphT(name, nodes, edges, lambda, numMarkers, numReads, 0);
		graph.setLocalRate(localRateMap);
    	return graph;
	}
	
	private static Alignments getConnection(Alignments left, Alignments right) {
		Alignments region=new Alignments(left.getChr(), left.getEnd(), right.getStart());
		region.setStrand("+");
		if(left.getEnd()>right.getEnd()){
			region=new Alignments(left.getChr(), right.getEnd(), left.getStart());
			region.setStrand("-");
		}
		return region;
	}

	private static EdgeSourceType getType(String line) {
		String[] tokens=line.split("\\[");
		if(tokens.length>1){
			String[] attributes=tokens[1].split(",");
			for(int i=0; i<attributes.length; i++){
				//System.err.println(attributes[i].split("=")[0]);
				if(attributes[i].split("=")[0].trim().equalsIgnoreCase("comment")){
					String type=attributes[i].split("=")[1];
					//System.err.println(type);
					if(type.equalsIgnoreCase(quote+"paired"+quote)){return EdgeSourceType.PAIRED;}
					if(type.equalsIgnoreCase("spliced")){return EdgeSourceType.SPLICED;}
				}
			}
		}
		return EdgeSourceType.SPLICED;
	}

	private static float getCount(String line) {
		String[] tokens=line.split("\\[");
		if(tokens.length>1){
			String[] attributes=tokens[1].split(",");
			for(int i=0; i<attributes.length; i++){
				if(attributes[i].split("=")[0].equals("label")){return new Float(attributes[i].split("=")[1]);}
			}
		}
		return 0;
	}
	
	private static float getNodeCount(String line) {
		String[] tokens=line.replaceAll("\\];", "").split("\\[");
		if(tokens.length>1){
			String[] attributes=tokens[1].split(",");
			for(int i=0; i<attributes.length; i++){
				if(attributes[i].split("=")[0].trim().equals("label")){
					return new Float(attributes[i].split("=")[1]);
				}
			}
		}
		return 0;
	}
	
	private static float getLocalRate(String line) {
		String[] tokens=line.replaceAll("\\];", "").split("\\[");
		if(tokens.length>1){
			String[] attributes=tokens[1].split(",");
			for(int i=0; i<attributes.length; i++){
				if(attributes[i].split("=")[0].trim().equals("comment")){
					return new Float(attributes[i].split("=")[1]);
				}
			}
		}
		return 0;
	}

	private static Alignments getNode(String node){
		node=node.split("\\[")[0];
		//System.err.println(node);
		Alignments nodeAlign=new Alignments(node.replaceAll(quote, ""));
		return nodeAlign;
	}
	
	private static String getStyle(String line){
		
		int startIndex=line.indexOf("[");
		int endIndex=line.indexOf("]");
		
		String attributesString=line.substring(startIndex+1, endIndex);
		//System.err.println(attributesString);
		
		String[] attributes=attributesString.split(",");
		for(int i=0; i<attributes.length; i++){
			if(attributes[i].split("=")[0].equals("style")){return attributes[i].split("=")[1];}
		}
		
		return "";
	}
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length==5 ){
			File dotFile=new File(args[0]);
			Alignments region=new Alignments(args[1], args[2], args[3]);
			String save=args[4];
			new ExtractRegionsFromDOT(dotFile, region, save);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=dot file \n args[1]=chr \n args[2]=start \n args[3]=end\n args[4] save file\n";
}
