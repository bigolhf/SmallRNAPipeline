/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import no.uio.medisin.bag.ngssmallrna.pipeline.IsomiRSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.MiRNAFeature;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 *   1. parse SAM file to extract and process the miRNA reads to determine isomiR content
 *   2. merge the counts from each sample to generate a single count file 
 *      that can be used for differential expression analysis
 * 
 *   Input is a SAM file
 * 
 * @author sr
 */

public class StepParseSAMForMiRNAs extends NGSStep{
    
    static Logger                       logger                      = LogManager.getLogger();
    static String                       FileSeparator               = System.getProperty("file.separator");
    
    private static final String         inFolder                    = "bowtie_genome_mapped";
    private static final String         miRNAAnalysisOutFolder      = "mirna_isomir_analysis";
    
    
    private static final String         infileExtension             = ".trim.clp.gen.sam";
    private static final String         isomirSummaryExtension      = ".trim.clp.gen.iso_summary.tsv";
    private static final String         isomirPrettyExtension       = ".trim.clp.gen.iso_pretty.tsv";
    private static final String         miRCountsExtension          = ".trim.clp.gen.mircounts.tsv";
    
    
    
    private StepInputData               stepInputData;
    private StepResultData              stepResultData;
    

    private List<MiRNAFeature>          miRBaseMiRNAList                   = new ArrayList<>();
    private List<MiRNAFeature>          miRNAHitList;
    private ArrayList<IsomiRSet>        isomiRList;
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepParseSAMForMiRNAs(StepInputData sid){
        stepInputData = sid;
    }
    
    @Override
    public void execute(){
        /*
            parseSAMmiRNAsParams.put("bleed", this.getSamParseForMiRNAsBleed());
            parseSAMmiRNAsParams.put("miRBaseHostGFFFile", this.getMiRBaseHostGFF());
            parseSAMmiRNAsParams.put("miRBaseRootFolder", this.getMirBaseVersionRoot());
            parseSAMmiRNAsParams.put("host", this.getBowtieMappingReferenceGenome());
            parseSAMmiRNAsParams.put("baseline_percent", this.getSamParseForMiRNAsBaselinePercent());
            parseSAMmiRNAsParams.put("analyze_isomirs", this.getSamParseForMiRNAsAnalyzeIsomirs());
        */
        try{
            stepInputData.verifyInputData();            
        }
        catch(IOException exIO){
            logger.info("exception parsing InputData" + exIO);
        }
    
        
        try{
            this.loadMiRBaseData((String) stepInputData.getStepParams().get("host"), (String) stepInputData.getStepParams().get("miRBaseHostGFFFile"));
        }
        catch(IOException ex){
            logger.error("error reading miRBase reference file <" + (String) stepInputData.getStepParams().get("miRBaseHostGFFFile") + ">\n" + ex.toString());
        }
        
        Boolean analyzeIsomirs = (Boolean) stepInputData.getStepParams().get("analyze_isomirs");

        String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
        String miRNAAnalysisOutputFolder = pathToData + FileSeparator + miRNAAnalysisOutFolder;
        miRNAAnalysisOutputFolder = miRNAAnalysisOutputFolder.replace(FileSeparator + FileSeparator, FileSeparator).trim();
        Boolean fA = new File(miRNAAnalysisOutputFolder).mkdir();       
        if (fA) logger.info("created output folder <" + miRNAAnalysisOutputFolder + "> for results" );
        String samLine = null;
        
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                
                int bleed = (int) stepInputData.getStepParams().get("bleed");
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                
                String samInputFile = pathToData + FileSeparator + inFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", infileExtension);
                logger.info(sampleData.getDataFile().replace(".fastq", infileExtension));
                int matchCount5 = 0;
                int matchCount3 = 0;
                int preMatchCount5 = 0;
                int preMatchCount3 = 0;
                int totalCounts = 0;
                samLine = null;
                BufferedReader brSAM = new BufferedReader(new FileReader(new File(samInputFile)));
                    isomiRList = new ArrayList<>();
                    miRNAHitList = new ArrayList<>();
                    while((samLine=brSAM.readLine())!= null){
                        /*
                            1   QNAME	   Query template NAME
                            2   FLAG	   bitwise FLAG
                            3   RNAME	   Reference sequence NAME
                            4   POS	   1-based leftmost mapping POSition
                            5   MAPQ	   MAPping Quality
                            6   CIGAR	   CIGAR string
                            7   RNEXT	   Ref. name of the mate/next read
                            8   PNEXT	   Position of the mate/next read
                            9   TLEN	   observed Template LENgth
                            10  SEQ	   segment SEQuence
                            11  QUAL	   ASCII of Phred-scaled base QUALity+33
                        
                        */
                        if(samLine.startsWith("@")) continue;
                        
                        totalCounts += Integer.parseInt(samLine.split("\t")[0].split("-")[1]);
                        if(samLine.split("\t")[1].equals("16") || samLine.split("\t")[1].equals("0")){
                            String strand = "";
                            if (samLine.split("\t")[1].equals("16")) {
                                strand = "-";
                                preMatchCount3++;
                            }
                            else{
                                strand = "+";
                                preMatchCount5++;
                            }
                            
                            int startPos = Integer.parseInt(samLine.split("\t")[3]);
                            String cigarStr = samLine.split("\t")[5].replace("M", "").trim();
                            int endPos = startPos + Integer.parseInt(cigarStr);
                            String chr = samLine.split("\t")[2].trim();
                            String mdString = samLine.split("\t")[12];
                            
                            MiRNAFeature miRNAFeature = this.doesReadOverlapKnownMiRNA(startPos, endPos, chr, strand, bleed);
                            if (miRNAFeature != null){
                                MiRNAFeature miRNAHit = new MiRNAFeature(miRNAFeature);
                                //logger.debug(miRNAHit.getName());
                                String name = samLine.split("\t")[0];
                                String sequence = samLine.split("\t")[9];
                                if(miRNAHitList.contains(miRNAHit)){ 
                                    miRNAHitList.get(miRNAHitList.indexOf(miRNAHit)).addIsomiR(name, startPos, cigarStr, mdString, sequence);
                                }
                                else{
                                    miRNAHit.addIsomiR(name, startPos, cigarStr, mdString, sequence);
                                    miRNAHitList.add(miRNAHit);
                                }
                                    
                                if(strand.equals("+")) matchCount5++;
                                else matchCount3++;
                                List<String> mutations = new ArrayList<>();
                                Matcher match = Pattern.compile("[0-9]+|[a-z]+|[A-Z]+").matcher(samLine.split("\t")[12].split(":")[2]);
                                String outputString = samLine.split("\t")[0] + ":" + startPos + ":" + endPos + ":[" + chr + "]:" + samLine.split("\t")[12] + ": ";
                                while (match.find()) {
                                    mutations.add(match.group());
                                    outputString = outputString.concat(match.group() + "|");
                                }
                                
                            }
                            
                            
                            
                            
                        }
                    }
                    logger.info("  total mapped counts = " + totalCounts);
                    Double minCounts = (double) totalCounts /100000.0;
                    logger.info((matchCount5 + matchCount3) + " reads (" + matchCount5 + " 5'" + "/" + matchCount3 + " 3' ) were mapped");
                    
                    if(analyzeIsomirs){
                        logger.info("  calculate isomiR dispersions");
                        for(MiRNAFeature miRHit: miRNAHitList){
                            if (miRHit.getTotalCounts() > minCounts.intValue()){
                                ArrayList isomirPtsAsHash = miRHit.characterizeIsomiRs((int) stepInputData.getStepParams().get("baseline_percent"), minCounts.intValue());
                                this.isomiRList.add(new IsomiRSet(miRHit.getMimatID(), sampleData.getNote(), sampleData.getDataFile().replace(".fastq", ""), isomirPtsAsHash));
                            }
                        }



                        logger.info("  write isomiRs");

                        String  isoDetailsFile = miRNAAnalysisOutputFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", isomirSummaryExtension);
                        String  isoPrettyFile  = miRNAAnalysisOutputFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", isomirPrettyExtension);

                        BufferedWriter brDetails = new BufferedWriter(new FileWriter(new File(isoDetailsFile)));
                        BufferedWriter brPretty  = new BufferedWriter(new FileWriter(new File(isoPrettyFile)));
                            for(MiRNAFeature miRHit: this.miRNAHitList){
                                if (miRHit.getTotalCounts() > minCounts.intValue()){
                                    logger.debug(miRHit.getName());
                                    brDetails.write(miRHit.reportIsomiRs((int) stepInputData.getStepParams().get("baseline_percent"), minCounts.intValue()));
                                    brPretty.write(miRHit.prettyReportIsomiRs((int) stepInputData.getStepParams().get("baseline_percent"), minCounts.intValue()));
                                }
                            }
                        brPretty.close();
                        brDetails.close();
                    }
                    
                    logger.info("  write miRNA counts");

                    String  miRCountsFile  = miRNAAnalysisOutputFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", miRCountsExtension);
                    
                    BufferedWriter brCounts  = new BufferedWriter(new FileWriter(new File(miRCountsFile)));
                        for(MiRNAFeature miR: this.miRBaseMiRNAList){
                            if(miRNAHitList.contains(miR)){
                                for(MiRNAFeature miRHit: this.miRNAHitList){
                                    if(miRHit.equals(miR)){
                                        brCounts.write(miR.getMimatID() + ":" + miR.getName() + "\t" + miRHit.getTotalCounts() + "\n");
                                        break;
                                    }
                                }
                            }
                            else{
                                brCounts.write(miR.getMimatID() + ":" + miR.getName() + "\t" + 0 + "\n");                                
                            }
                        }
                    brCounts.close();
                    
                brSAM.close();
                logger.info("  completed processing SAM file\n\n");
                
                
            }
            catch(IOException ex){
                logger.error("error executing parse SAM for miRNAs command\n" + ex.toString());
            }
            catch(ArrayIndexOutOfBoundsException exBnd){
                logger.error("error parsing line " + samLine);
                logger.error(exBnd);
            }
        }
        
        if(analyzeIsomirs){
            String dispersionFile   = miRNAAnalysisOutputFolder + FileSeparator + stepInputData.getProjectID() + ".disp.tsv";
            String summaryFile      = miRNAAnalysisOutputFolder + FileSeparator + stepInputData.getProjectID() + ".disp.summary.tsv";
            logger.info("write dispersions to file <" + dispersionFile + ">");
            try{
                BufferedWriter bwDp = new BufferedWriter(new FileWriter(new File(dispersionFile)));
                BufferedWriter bwSm = new BufferedWriter(new FileWriter(new File(summaryFile)));   
                    bwSm.write(IsomiRSet.printSummaryHeader());
                    for(IsomiRSet isomiRset: isomiRList){
                        isomiRset.calcDistParameters();
                        bwSm.write(isomiRset.printSummary());
                        bwDp.write(isomiRset.tabReportIsomiRSet());                
                    }
                bwSm.close();
                bwDp.close();
            }
            catch(IOException exIO){
                logger.info("error writing isomiR dispersion File <" + dispersionFile + ">\n" + exIO);
            }
        }
        
        
    }
    
    
    
    
    /**
     * Does the read sufficiently overlap a defined miRNA entry?
     * 
     * @param start
     * @param stop
     * @param chr
     * @param strand
     * @param bleed     int : specifies how much a read can 'miss' an entry
     *                        and still be counted
     * 
     * @return MiRNAFeature
     */
    public MiRNAFeature doesReadOverlapKnownMiRNA(int start, int stop, String chr, String strand, int bleed){
        
        for(MiRNAFeature miRBaseEntry: this.miRBaseMiRNAList){
            if (miRBaseEntry.chromosomeMatch(chr)){
                if(strand.equals(miRBaseEntry.getStrand())){
                    
                    if( java.lang.Math.abs(start - miRBaseEntry.getStartPos()) <= bleed){

                        if( java.lang.Math.abs(stop - miRBaseEntry.getEndPos()) <= bleed){
                            return miRBaseEntry;
                        }

                    }

                }
            }
        }
        
        return null;
        
    }
    
    
    
    
    /**
     * Verify Input Data for parsing SAM file for miRNAs
     * 
     */        
    @Override
    public void verifyInputData(){
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            /*
            if (sampleData.getDataFile().toUpperCase().endsWith(infileExtension.toUpperCase())==false)
            {
                throw new IllegalArgumentException("AdapterTrimming: incorrect file extension for input file <" 
                  + sampleData.getDataFile() + ">. " 
                  + "should have <" + infileExtension + "> as extension");
            }
            
            if (sampleData.getDataFile().toUpperCase().endsWith(outfileExtension.toUpperCase())==true)
            {
                logger.warn("AdapterTrimming: input file has output file extension (.trim.fastq)");
                logger.warn("this file has already been trimmed");
            }
            */
            
        }
            // does input file have correct extension?
        // does input file have the same extension as expected for the output file?
    }
    
    /**
     * 1. 
     * load miRNA specs (name and Chromosome position) from GFF file
     * downloaded from miRBase
     * Because different releases of miRBase use different releases of
     * reference genome, we have to track both miRBase and genome reference IDs
     * 2.
     * Load Sequence from Mature.fa file
     * 
     * @param host              String : 3 char abbreviation for host
     * @param miRBaseGFFFile    String : absolute path to file
     * @throws IOException
     * 
     */
    public void loadMiRBaseData(String host, String miRBaseGFFFile) throws IOException{

        HashMap <String, String> miRBaseSeq = new HashMap();
        String matureFAFile = miRBaseGFFFile.replace("gff3", "mature.fa");
        //String matureFAFile = new File(miRBaseGFFFile).getParent() + FileSeparator + "mature.fa";
        BufferedReader brFA = new BufferedReader(new FileReader(new File(matureFAFile)));
        String lineFA = null;
        while ((lineFA = brFA.readLine())!=null){

            String seq = brFA.readLine().trim();
            String entryHost = lineFA.split(" ")[0].substring(1).split("-")[0].trim();
            if(entryHost.equals(host)){
                String mimatID = lineFA.split(" ")[1].trim();
                miRBaseSeq.put(mimatID, seq);
            }
            
        }
        
        
        
        String line = null;
        BufferedReader brMiR = new BufferedReader(new FileReader(new File(miRBaseGFFFile)));
            while((line = brMiR.readLine())!= null){
                
                if(line.startsWith("#")) continue;
                if(line.contains("miRNA_primary_transcript")) continue;
                /*
                    chr1            chromosome
                    source          n/a here
                    miRNA           feature type (n/a)
                    start pos
                    end pos
                    score           n/a here               
                    strand          (+/-)
                    frame           n/a here
                    attributes      e.g. ID=MIMAT0027619;Alias=MIMAT0027619;Name=hsa-miR-6859-3p;Derives_from=MI0022705
                
                */
                String chr = line.split("\t")[0].trim();
                if(chr.contains("chr")) chr = chr.replace("chr", "");
                int startPos = Integer.parseInt(line.split("\t")[3].trim());
                int endPos = Integer.parseInt(line.split("\t")[4].trim());
                String strand = line.split("\t")[6].trim();
                
                String id = "";
                String name = "";
                String parent = "";
                String attribs[] = line.split("\t")[8].split(";");
                
                for (String attribStr: attribs){
                    String attribType = attribStr.split("=")[0].trim();
                    String attribValue = attribStr.split("=")[1].trim();
                    switch (attribType){
                        case "ID":
                            id = attribValue;
                            break;
                            
                        case "Alias":                            
                            break;
                            
                        case "Name":
                            name = attribValue;
                            break;
                            
                        case "Derives_from":
                            parent = attribValue;
                            break;
                            
                        default:
                            logger.warn("unknown attribute in parsing miRNA entry from GFF file " + miRBaseGFFFile + ">");
                            break;
                    }
                }
                // this may need revising. In some cases we dont care if there is no sequence, we still want the feature
                String seq = miRBaseSeq.get(id);
                this.miRBaseMiRNAList.add(new MiRNAFeature(name, chr, startPos, endPos, strand, id, parent, seq));
                if(seq == null) 
                    logger.warn("no sequence found for entry <" + id + ">. Skipping");

            }
        brMiR.close();
        logger.info("read " + miRBaseMiRNAList.size() + " miRNA entries");
        
    }
    
    
    
    @Override
    public void outputResultData(){
        
    }
}
