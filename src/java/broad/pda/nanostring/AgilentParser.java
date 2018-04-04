package broad.pda.nanostring;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.pda.differentialExpression.DifferentialExpression;
import broad.pda.geneexpression.ExpressionExperimentInfo;
import broad.pda.geneexpression.agilent.AgilentUtils;

public class AgilentParser {
	
	Collection<String> positives;
	Collection<String> negatives;
	double floorVal=0;
	private double foldCutoff=1.0;
	private boolean useFold=false;
	private int numPerm=1000;
	private double[] fudgeFactors={0,.01,.1,1};
	private double alpha=0.05;

	//Parse file and normalize
	public AgilentParser(File file, File description, String save) throws IOException, ParseException{
		MatrixWithHeaders data=parseNanostring(file);
		
		Map<String, Collection<String>> experimentInfo=AgilentUtils.parseExperimentInfoFileToGroups(description);
		Map<String, ExpressionExperimentInfo> fullInfo=AgilentUtils.parseExperimentInfoFile(description);
		
		//normalizeByMedianPositiveSpikeIns(data);
		//normalizeByExperimentMedian(data);
		//floorNonExpressed(data);
		
		//filter non-expressed genes
		//data=removeNonExpressed(data);
		
		//Convert to log2 scale
		//convertToLog(data);
		
		//Median normalize
		//medianNorm(data);
		
		//Compute Z-Scores within batch and compute significance off of these
		computeZScoresByBatch(data, fullInfo, experimentInfo);
		
		
		//compute correlation of replicates
		for(String group: experimentInfo.keySet()){
			Collection<String> samples=experimentInfo.get(group);
			MatrixWithHeaders subdata=data.submatrixByColumnNames(samples);
			for(int i=0; i<subdata.columnDimension(); i++){
				for(int j=i+1; j<subdata.columnDimension(); j++){
					double r=Statistics.pearsonDistance(subdata.getColumn(i), subdata.getColumn(j));
					System.out.println(group+" "+subdata.getColoumnName(i)+" "+subdata.getColoumnName(j)+" "+r);
				}
			}
		}
		
		
		//Compute fold change
		MatrixWithHeaders foldMatrix=foldChange(data, experimentInfo);
		foldMatrix.writeGCT(save+".fold.gct");
		
		
		//TODO discretize matrix based on fold cutoff (average, and min fold)
		
		//TODO Compute FDR values using permutations, discretize using FDR
		MatrixWithHeaders fdrMatrix=computeFDRMatrix(data, experimentInfo);
		fdrMatrix.writeGCT(save+".fdr.gct");
		
		data.writeGCT(save, experimentInfo, fullInfo, true);
	}

	
	private MatrixWithHeaders computeFDRMatrix(MatrixWithHeaders data,	Map<String, Collection<String>> experimentInfo) {
		ArrayList<String> includedGroups=new ArrayList(experimentInfo.keySet());
		includedGroups.remove("Control");
		MatrixWithHeaders rtrn=new MatrixWithHeaders(data.getRowNames(), includedGroups);
		Collection<String> controls=experimentInfo.get("Control");
		
		for(String group: experimentInfo.keySet()){
			if(!group.equalsIgnoreCase("Control")){
				Collection<String> group1=experimentInfo.get(group);
				DifferentialExpression diff=new DifferentialExpression(data, group1, controls, useFold, false, numPerm, fudgeFactors, null, null);
				MatrixWithHeaders fdr=diff.getFDRMatrix();//by sign
				MatrixWithHeaders fdrAbs=diff.getAbsFDRMatrix(); //by absolute val
				for(String gene: data.getRowNames()){
					//double fdrVal=Math.min(Statistics.min(fdr.getRow(gene)), Statistics.min(fdrAbs.getRow(gene)));
					double fdrVal=Statistics.min(fdr.getRow(gene));
					rtrn.set(gene, group, fdrVal);
				}
			}
		}
		
		
		
		return rtrn;
	}

	private void computeZScoresByBatch(MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> fullInfo, Map<String, Collection<String>> experimentInfo) throws IOException {
		//Collection<String> controls=experimentInfo.get("Control");
		//Map<String, Collection<String>> controlsByBatch=getControlsByBatch(controls, fullInfo);
		
		/*for(String group: experimentInfo.keySet()){
			Collection<String> experiments=experimentInfo.get(group);
			for(String experiment: experiments){
				ExpressionExperimentInfo info=fullInfo.get(experiment);
				Collection<String> batchControls=controlsByBatch.get(info.getBatch());
				MatrixWithHeaders batchControlData=data.submatrixByColumnNames(batchControls);
				for(String gene: data.getRowNames()){
					double val=data.get(gene, experiment);
					double zScore=zScore(val, batchControlData.getRow(gene));
					data.set(gene, experiment, zScore);
				}
			}  
		}*/
		
		Map<String, Collection<String>> experimentsByBatch=new TreeMap<String, Collection<String>>();
		for(String group: experimentInfo.keySet()){
			Collection<String> experiments=experimentInfo.get(group);
			for(String experiment: experiments){
				String batch=fullInfo.get(experiment).getBatch();
				Collection<String> all=new TreeSet<String>();
				if(experimentsByBatch.containsKey(batch)){all=experimentsByBatch.get(batch);}
				all.add(experiment);
				experimentsByBatch.put(batch, all);
			}
		}
		
		for(String batch: experimentsByBatch.keySet()){
			MatrixWithHeaders batchMatrix=data.submatrixByColumnNames(experimentsByBatch.get(batch));
			
			Collection<String> controls=new TreeSet<String>();
			for(String experiment: batchMatrix.getColumnNames()){
				if(fullInfo.get(experiment).isControl()){controls.add(experiment);}
			}
			
			for(String experiment: batchMatrix.getColumnNames()){
				for(String gene: batchMatrix.getRowNames()){
					double val=batchMatrix.get(gene, experiment);
					double zScore=zScore(val, data.getValues(gene, controls));
					data.set(gene, experiment, zScore);
				}
			}
			
			batchMatrix.writeGCT("batch"+batch+".gct");
		}
	}

	private Map<String, Collection<String>> getControlsByBatch(Collection<String> controls,	Map<String, ExpressionExperimentInfo> fullInfo) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String experiment: controls){
			String batch=fullInfo.get(experiment).getBatch();
			Collection<String> all=new TreeSet<String>();
			if(rtrn.containsKey(batch)){all=rtrn.get(batch);}
			all.add(experiment);
			rtrn.put(batch, all);
		}
		
		
		return rtrn;
	}

	private double zScore(double val, double[] row) {
		return Statistics.zScore(val, row, .1);
	}

	private MatrixWithHeaders foldChange(MatrixWithHeaders data, Map<String, Collection<String>> experimentInfo) {
		ArrayList<String> experiments=getExperiments(experimentInfo);
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(data.getRowNames(), experiments);
		
		Collection<String> controls=experimentInfo.get("Control");
		
		MatrixWithHeaders controlData=data.submatrixByColumnNames(controls);
		
		Collection<String> foldList=new TreeSet<String>();
		
		for(String group: experimentInfo.keySet()){
				//if(!group.equalsIgnoreCase("Control")){
				System.err.println(group);
				for(String experiment: experimentInfo.get(group)){
					//MatrixWithHeaders sampleData=data.submatrixByColumnNames(experimentInfo.get(group));
					for(String gene: data.getRowNames()){
						double fold=fold(data.get(gene, experiment), controlData.getRow(gene));
						rtrn.set(gene, experiment, fold);
					}
					//if(fold>foldCutoff){foldList.add(gene);}
				}
			//}
		}
		
		System.err.println("Filtered fold "+foldList.size()+" "+data.getRowNames().size());
		return rtrn;
	}
	
	
	private double fold(double val, double[] values) {
		return val-Statistics.average(values);
	}
	private ArrayList<String> getExperiments(Map<String, Collection<String>> experimentInfo) {
		ArrayList<String> rtrn=new ArrayList<String>();
		
		for(String group: experimentInfo.keySet()){
			Collection<String> experiments=experimentInfo.get(group);
			rtrn.addAll(experiments);
		}
		
		return rtrn;
	}
	
	private MatrixWithHeaders computeTStats(MatrixWithHeaders data, Map<String, Collection<String>> experimentInfo){
		MatrixWithHeaders rtrn=new MatrixWithHeaders(data.getRowNames(), new ArrayList(experimentInfo.keySet()));
		
		Collection<String> controls=experimentInfo.get("Control");
		
		Collection<String> foldList=new TreeSet<String>();
		
		for(String group: experimentInfo.keySet()){
			//if(!group.equalsIgnoreCase("Control")){
				System.err.println(group);
				for(String gene: data.getRowNames()){
					double tstat=tstat(data, gene, controls, experimentInfo.get(group));
					rtrn.set(gene, group, tstat);
				}
			//}
		}
		
		//System.err.println("Filtered fold "+foldList.size()+" "+data.getRowNames().size());
		return rtrn;
	}


	private double fold(double[] controlVals, double[] sampleVals) {
		return Statistics.median(sampleVals)-Statistics.median(controlVals);
	}
	
	private double fold(MatrixWithHeaders data, String gene, Collection<String> controls, String sample) {
		double[] controlVals=data.submatrixByColumnNames(controls).getRow(gene);
		double sampleVal=data.get(gene, sample);
		
		return sampleVal-Statistics.median(controlVals);
	}
	
	private double tstat(MatrixWithHeaders data, String gene, Collection<String> controls, Collection<String> samples) {
		double[] controlVals=data.submatrixByColumnNames(controls).getRow(gene);
		double[] sampleVals=data.submatrixByColumnNames(samples).getRow(gene);
		
		if(sampleVals.length==1){return Statistics.zScore(sampleVals[0], controlVals,.1);}
		return Statistics.tstat(sampleVals, controlVals,.1);
	}


	private MatrixWithHeaders removeNonExpressed(MatrixWithHeaders data) {
		Collection<String> expressed=new TreeSet<String>();
		
		for(String gene: data.getRowNames()){
			double[] vals=data.getRow(gene);
			int counter=0;
			for(int i=0; i<vals.length; i++){
				if(vals[i]>2){counter++;}
			}
			if(counter>1){expressed.add(gene);}
		}
		
		System.err.println("Expressed "+expressed.size()+" total "+data.getRowNames().size());
		return data.submatrixByRowNames(expressed);
	}


	private void medianNorm(MatrixWithHeaders data) {
		for(String experiment: data.getColumnNames()){
			double median=Statistics.median(data.getColumn(experiment));
			for(String gene: data.getRowNames()){
				double val=data.get(gene, experiment);
				double norm=val-median;
				data.set(gene, experiment, norm);
			}
		}
		
	}


	private void convertToLog(MatrixWithHeaders data) {
		for(String gene: data.getRowNames()){
			for(String experiment: data.getColumnNames()){
				double val=data.get(gene, experiment);
				if(val==0){System.err.println("HAS 0 value"); val=.1;}
				double log=Math.log(val)/Math.log(2);
				data.set(gene, experiment, log);
			}
		}
		
	}


	private void floorNonExpressed(MatrixWithHeaders data) {
		MatrixWithHeaders negativeProbes=getNegativeProbes(data);
		
		//Actually, I'm going to give it a fold over the max negatives. This way all values are non-zero and positive and easy to interpret
		for(String experiment: data.getColumnNames()){
			double[] negatives=negativeProbes.getColumn(experiment);
			double max=Statistics.max(negatives);
			for(String gene: data.getRowNames()){
				double val=data.get(gene, experiment);
				double norm=Math.max(val/max,1);
				data.set(gene, experiment, norm);
			}
		}
		
	}


	private void normalizeByExperimentMedian(MatrixWithHeaders data) {
		MatrixWithHeaders sampleProbes=getSampleProbes(data);
		
		for(String experiment: data.getColumnNames()){
			double[] vals=sampleProbes.getColumn(experiment);
			double median=Statistics.median(vals);
			
			for(String gene: data.getRowNames()){
				double val=data.get(gene, experiment);
				double norm=val/median;
				data.set(gene, experiment, norm);
			}
			
		}
		
	}

	
	private MatrixWithHeaders getNegativeProbes(MatrixWithHeaders data) {
		Collection<String> positives=new ArrayList<String>();
		
		Map<String, String> probeClass=data.getNanostringProbeClasses();
		for(String gene: data.getRowNames()){
			if(probeClass.get(gene).equalsIgnoreCase("Negative")){positives.add(gene);}
		}
		
		return data.submatrixByRowNames(positives);
	}

	private MatrixWithHeaders getSampleProbes(MatrixWithHeaders data) {
		Collection<String> positives=new ArrayList<String>();
		
		Map<String, String> probeClass=data.getNanostringProbeClasses();
		for(String gene: data.getRowNames()){
			if(probeClass.get(gene).equalsIgnoreCase("Endogenous")){positives.add(gene);}
		}
		
		return data.submatrixByRowNames(positives);
	}


	private void normalizeByMedianPositiveSpikeIns(MatrixWithHeaders data) {
		
		MatrixWithHeaders positives=getPositives(data);
		
		Collection<Double> vals=new ArrayList<Double>();
		for(String experiment: positives.getColumnNames()){
			double sum=Statistics.sum(positives.getColumn(experiment));
			vals.add(sum);
		}
		
		double median=Statistics.median(vals);
		
		//Normalize sample values
		for(String experiment: data.getColumnNames()){
			double sum=Statistics.sum(positives.getColumn(experiment));
			double experimentNormFactor=median/sum;
			for(String gene: data.getRowNames()){
				double val=data.get(gene, experiment);
				double norm=val*experimentNormFactor;
				data.set(gene, experiment, norm);
			}
		}
		
	}


	private MatrixWithHeaders getPositives(MatrixWithHeaders data) {
		Collection<String> positives=new ArrayList<String>();
		
		Map<String, String> probeClass=data.getNanostringProbeClasses();
		for(String gene: data.getRowNames()){
			if(probeClass.get(gene).equalsIgnoreCase("Positive")){positives.add(gene);}
		}
		
		return data.submatrixByRowNames(positives);
	}


	private MatrixWithHeaders parseNanostring(File file) throws IOException, ParseException{
		return new MatrixWithHeaders(file.getAbsolutePath());
	}
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>1){
			File file=new File(args[0]);
			File description=new File(args[1]);
			String save=args[2];
			new AgilentParser(file, description, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=file \n args[1]=description \n args[2]=save";
	
}
