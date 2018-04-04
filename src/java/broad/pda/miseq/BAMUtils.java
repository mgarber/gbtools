package broad.pda.miseq;

import java.io.IOException;

import broad.core.util.PipelineUtils;

public class BAMUtils {

	public static void markDuplicates(String bamFile, String save, String info, String bsub, Runtime run, String jobID) throws IOException, InterruptedException{
		String cmd="java -jar /seq/mgarber/tools/picard-tools-1.66/MarkDuplicates.jar";
		cmd+=" I="+bamFile+" O="+save+" M="+info;
		PipelineUtils.bsubProcess(run, jobID, cmd, bsub);
	}
	
	public static void buildBAMIndex(){
		
	}
	
}
