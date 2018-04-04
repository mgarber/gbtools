package broad.core.alignment;

import java.util.HashMap;
import java.util.Iterator;

public class AlignmentStats {
	private HashMap<String,long[]> queryStats   = new HashMap<String,long[]>(); 
	private HashMap<String,long[]> subjectStats = new HashMap<String,long[]>();
	long maxQuery;
	long maxSubject;
	long minQuery;
	long minSubject;
	
	public void add(AlignmentSummary summary) {
		String subject = summary.getSubject();
		String query   = summary.getQuery();
		
		long[] sStats  = subjectStats.get(subject);
		long[] qStats  = queryStats.get(query);
		//System.out.println(summary.toString());
		maxQuery = Math.max(maxQuery,summary.getQueryEnd());
		minQuery = minQuery == 0 ? summary.getQueryStart() : Math.min(minQuery,summary.getQueryStart());
		maxSubject = Math.max(maxSubject, summary.getSubjectEnd());
		minSubject = minSubject == 0 ? summary.getSubjectStart() : Math.min(minSubject,summary.getSubjectStart());
		
		if(sStats == null){
			sStats = new long[2];
			sStats[0] = summary.getSubjectStart();
			sStats[1] = summary.getSubjectEnd();
			subjectStats.put(subject,sStats);
		} else {
			sStats[0] = Math.min(sStats[0],summary.getSubjectStart());
			sStats[1] = Math.max(sStats[1],summary.getSubjectEnd());
		}
		if(qStats == null){
			qStats = new long[2];
			qStats[0] = summary.getQueryStart();
			qStats[1] = summary.getQueryEnd();
			queryStats.put(query,qStats);
		} else {
			qStats[0] = Math.min(qStats[0],summary.getQueryStart());
			qStats[1] = Math.max(qStats[1],summary.getQueryEnd());
		}
		
	}
	
	public HashMap<String,long[]> getSubjectStats() {
		return subjectStats;
	}
	
	public HashMap<String,long[]> getQueryStats() {
		return queryStats;
	}
	
	public long getMinQuery() {
		return minQuery;			
	}
	
	public long getMaxQuery() {
		return maxQuery;
	}
	
	public long getMinSubject() {
		return minSubject;
	}
	
	public long getMaxSubject() {
		return maxSubject;
	}
	
	public String toString() {
		StringBuffer buf = new StringBuffer("Alignment Summary Statistics");
		buf.append("maxSubj<").append(maxSubject)
			.append("> minSubject<").append(minSubject)
			.append("> maxQuery<").append(maxQuery)
			.append("> minQuery<").append(minQuery).append(">\n\t");
			
		Iterator<String> it = queryStats.keySet().iterator();
		buf.append("Query:\n");
		while(it.hasNext()) {
			String query = it.next();
			buf.append("\t\tName ").append(query);
			buf.append(" min ").append(queryStats.get(query)[0])
				.append(" max ").append(queryStats.get(query)[1])
				.append("\n");
		}
		
		it = subjectStats.keySet().iterator();
		buf.append("\tSubject:\n");
		while(it.hasNext()) {
			String subject = it.next();
			buf.append("\t\tName ").append(subject);
			buf.append(" min ").append(subjectStats.get(subject)[0])
				.append(" max ").append(subjectStats.get(subject)[1])
				.append("\n");
		}
		
		return buf.toString();
		
	}
	
}
