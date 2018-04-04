package broad.core.alignment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import broad.core.annotation.LightweightGenomicAnnotation;




public class SmatchAlignmentReader extends AlignmentList<SmatchAlignmentSummary> {

    private static final String MATCHING = "MATCHING";

    private static final String SEQUENCE = "\tSequence:";

    public SmatchAlignmentReader(String alignmentFile) throws IOException {
        super();
        loadAll(alignmentFile);
    }

    public SmatchAlignmentReader() {
        super();
    }

    public void loadAll(String alignmentFile) throws IOException {
        loadAll(alignmentFile, 0);
    }

    public void loadAll(String alignmentFile, int shiftQuery) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(alignmentFile));

        String line;
        String[] primerData = null;
        while((line = br.readLine()) != null) {
            if(line.startsWith("#") || line.trim().length() ==0)
                continue;

            if (line.startsWith(MATCHING)) {

                primerData = line.split(" +");
                if (primerData.length < 3)
                    throw new IOException("The line " + line + " has fewer than 3 tokens");

            } else if (line.startsWith(SEQUENCE) && primerData != null) {

                String [] queryData = line.split(" +");
                if (queryData.length < 12)
                    throw new IOException("The line " + line + " has fewer than 12 tokens");

                SmatchAlignmentSummary summary = new SmatchAlignmentSummary(primerData, queryData);

                LightweightGenomicAnnotation query = summary.getA();
                query.setStart(query.getStart() + shiftQuery);
                query.setEnd(query.getEnd() + shiftQuery);

                addAlignment(summary);
            }
        }

        br.close();
    }


}
