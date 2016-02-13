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
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import no.uio.medisin.bag.ngssmallrna.pipeline.IsomiRSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.MiRNAFeature;
import no.uio.medisin.bag.ngssmallrna.pipeline.MirFeatureSet;
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

public class StepParseSAMForMiRNAs extends NGSStep implements NGSBase{
    
    static Logger                   logger                          = LogManager.getLogger();
    
    public  static final String     STEP_ID_STRING                  = "ParseSAMForMiRNAs";
    private static final String     ID_BLEED                        = "bleed";
    private static final String     ID_ISOMIRS                      = "analyzeIsomirs";
    private static final String     ID_MIRBASE_VERSION              = "mirbaseVersion";
    private static final String     ID_REF_GENOME                   = "host";
    private static final String     ID_BASELINE                     = "baselinePercent";
        

    private static final String     INFILE_EXTENSION                = ".trim.clp.gen.sam";
    private static final String     ISOMIR_SUMMARY_EXTENSION        = ".trim.clp.gen.iso_summary.tsv";
    private static final String     ISOMIR_PRETTY_EXTENSION         = ".trim.clp.gen.iso_pretty.tsv";
    private static final String     MIRCOUNTS_EXTENSION             = ".trim.clp.gen.mircounts.tsv";

    private List<MiRNAFeature>      miRNAHitList;
    private ArrayList<IsomiRSet>    isomiRList;
    MirFeatureSet                   mirBaseSet                      = new MirFeatureSet();           
    
    private int                     locationBleed                   = 2;
    private Boolean                 analyzeIsomirs                  = false;
    private int                     miRBaseRelease                  = 20;
    private String                  referenceGenome                 = "";
    private int                     baselinePercent                 = 5;
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepParseSAMForMiRNAs(StepInputData sid){
        stepInputData = sid;
    }
    
    
    
    
    
    /**
     * in this method we are simply checking that the configuration file 
     * has all the entries we need. We dont check if the values are acceptable
     * that is the role of the NGSStep.
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{
        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        if(configData.get(ID_BLEED)==null) {
            throw new NullPointerException("<" + ID_BLEED + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_ISOMIRS)==null) {
            throw new NullPointerException("<" + ID_ISOMIRS + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIRBASE_VERSION)==null) {
            throw new NullPointerException("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_BASELINE)==null) {
            throw new NullPointerException("<" + ID_BASELINE + "> : Missing Definition in Configuration File");
        }
        

      
        try{
            this.setMiRBaseRelease((Integer) configData.get(ID_MIRBASE_VERSION));
        }
        catch(Exception exNm){
            throw new NumberFormatException(ID_MIRBASE_VERSION + " <" + configData.get(ID_MIRBASE_VERSION) + "> is not an integer");
        }        
        if (this.getMiRBaseRelease() <= 0){
            throw new IllegalArgumentException(ID_MIRBASE_VERSION + " <" + configData.get(ID_MIRBASE_VERSION) + "> must be positive integer");
        }

        
        try{
            this.setBaselinePercent((Integer) configData.get(ID_BASELINE));
        }
        catch(Exception exNm){
            throw new NumberFormatException(ID_BASELINE + " <" + configData.get(ID_BASELINE) + "> is not an integer");
        }        
        if (this.getBaselinePercent() <= 0){
            throw new IllegalArgumentException(ID_BASELINE + " <" + configData.get(ID_BASELINE) + "> must be an integer between 0 and 100");
        }

        try{
            this.setLocationBleed((Integer) configData.get(ID_BLEED));
        }
        catch(Exception exNm){
            throw new NumberFormatException(ID_BLEED + " <" + ID_BLEED + "> is not an integer");
        }        
        if (this.getLocationBleed() <= 0 ){
            throw new IllegalArgumentException(ID_BLEED + " <" + configData.get(ID_BLEED) + "> must be > 0 ");
        }        

        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }

        try{
            this.setAnalyzeIsomirs((Boolean) configData.get(ID_ISOMIRS));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_BLEED + " <" + configData.get(ID_BLEED) + "> cannot be cast as Boolean");
        }        
        

        logger.info("passed");
    }
    
    
    
    /**
     * count up reads that overlap features specified in the GFF file
     * 
     * @throws IOException 
     */
    @Override
    public void execute()  throws IOException{
        
        
        this.setPaths();
        stepInputData.verifyInputData();            
    
        String gffFileMirBase = stepInputData.getDataLocations().getMirbaseFolder() + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3";
        String faFileMirBase = gffFileMirBase.replace("gff3", "mature.fa");
        mirBaseSet.loadMiRBaseData(this.getReferenceGenome(),gffFileMirBase, faFileMirBase);
        
        Boolean fA = new File(outFolder).mkdir();       
        if (fA) logger.info("created output folder <" + outFolder + "> for results" );
        String samLine = null;
        String samInputFile = "";
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                
                int bleed = this.getLocationBleed();
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                
                samInputFile = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION);
                logger.info(sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION));
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
                    /*
                        the following is rather approximate.
                        apparently, for 5,000,000 reads, the lowest detectable by qPCR is 50. so, we divide total counts by 100000
                        there has to be a better way....
                    */
                    Double minCountsForSingleFeature = (double) totalCounts /100000.0; // <= this is rather approximate
                    logger.info((matchCount5 + matchCount3) + " reads (" + matchCount5 + " 5'" + "/" + matchCount3 + " 3' ) were mapped");
                    
                    if(analyzeIsomirs){
                        logger.info("  calculate isomiR dispersions");
                        for(MiRNAFeature miRHit: miRNAHitList){
                            if (miRHit.getTotalCounts() > minCountsForSingleFeature.intValue()){
                                ArrayList isomirPtsAsHash = miRHit.characterizeIsomiRs(this.getBaselinePercent());
                                this.isomiRList.add(new IsomiRSet(miRHit.getMimatID(), sampleData.getNote(), sampleData.getFastqFile1().replace(".fastq", ""), isomirPtsAsHash));
                            }
                        }



                        logger.info("  write isomiRs");

                        String  isoDetailsFile = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", ISOMIR_SUMMARY_EXTENSION);
                        String  isoPrettyFile  = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", ISOMIR_PRETTY_EXTENSION);

                        BufferedWriter brDetails = new BufferedWriter(new FileWriter(new File(isoDetailsFile)));
                        BufferedWriter brPretty  = new BufferedWriter(new FileWriter(new File(isoPrettyFile)));
                            for(MiRNAFeature miRHit: this.miRNAHitList){
                                if (miRHit.getTotalCounts() > minCountsForSingleFeature.intValue()){
                                    logger.debug(miRHit.getName());
                                    brDetails.write(miRHit.reportIsomiRs(this.getBaselinePercent(), minCountsForSingleFeature.intValue()));
                                    brPretty.write(miRHit.prettyReportIsomiRs(this.getBaselinePercent(), minCountsForSingleFeature.intValue()));
                                }
                            }
                        brPretty.close();
                        brDetails.close();
                    }
                    
                    logger.info("  write miRNA counts");

                    String  miRCountsFile  = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", MIRCOUNTS_EXTENSION);
                    
                    BufferedWriter brCounts  = new BufferedWriter(new FileWriter(new File(miRCountsFile)));
                        for(MiRNAFeature miR: this.mirBaseSet.getMiRBaseMiRNAList()){
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
                logger.error("error processing sample <" + samInputFile + ">\n" + ex.toString());
                throw new IOException(STEP_ID_STRING + ": error processing sample <" + samInputFile + ">");
            }
            catch(ArrayIndexOutOfBoundsException exBnd){
                logger.error("error parsing line " + samLine);
                logger.error(exBnd);
                throw new IOException(STEP_ID_STRING + ": error processing sample <" + samInputFile + ">: samLine was \n" + samLine);
            }
        }
        
        if(analyzeIsomirs){
            String dispersionFile   = outFolder + FILESEPARATOR + stepInputData.getProjectID() + ".disp.tsv";
            String summaryFile      = outFolder + FILESEPARATOR + stepInputData.getProjectID() + ".disp.summary.tsv";
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
                throw new IOException(STEP_ID_STRING + "error writing isomiR dispersion File <" + dispersionFile + ">");
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
        
        for(MiRNAFeature miRBaseEntry: this.mirBaseSet.getMiRBaseMiRNAList()){
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
    public void verifyInputData() throws IOException{
        
        this.setPaths();
        
        // check the SAM files exist
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            if (sampleData.getFastqFile1()==null) throw new IOException("no Fastq1 file specified");
            String samInputFile = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION);            
            
            if ((new File(samInputFile)).exists()==false){
                throw new IOException(STEP_ID_STRING + ": SAM file <" + samInputFile + "> does not exist");
            }
            if (samInputFile.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + samInputFile + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
        }

    }
    
    
    
    
    /**
     * generate sample configuration data so the user can see what can be
     * specified
     *
     * @return
     */
    @Override
    public HashMap generateExampleConfigurationData() {

        logger.info(STEP_ID_STRING + ": generate example configuration data");

        HashMap configData = new HashMap();

        configData.put(ID_REF_GENOME, "hsa");
        configData.put(ID_BLEED, 2);
        configData.put(ID_BASELINE, 5);
        configData.put(ID_MIRBASE_VERSION, 20);
        configData.put(ID_ISOMIRS, true);

        return configData;
        
    }





    @Override
    public void verifyOutputData(){
        
    }

    /**
     * @return the locationBleed
     */
    public int getLocationBleed() {
        return locationBleed;
    }

    /**
     * @param locationBleed the locationBleed to set
     */
    public void setLocationBleed(int locationBleed) {
        this.locationBleed = locationBleed;
    }

    /**
     * @return the analyzeIsomirs
     */
    public Boolean getAnalyzeIsomirs() {
        return analyzeIsomirs;
    }

    /**
     * @param analyzeIsomirs the analyzeIsomirs to set
     */
    public void setAnalyzeIsomirs(Boolean analyzeIsomirs) {
        this.analyzeIsomirs = analyzeIsomirs;
    }

    /**
     * @return the miRBaseRelease
     */
    public int getMiRBaseRelease() {
        return miRBaseRelease;
    }

    /**
     * @param miRBaseRelease the miRBaseRelease to set
     */
    public void setMiRBaseRelease(int miRBaseRelease) {
        this.miRBaseRelease = miRBaseRelease;
    }

    /**
     * @return the ReferenceGenome
     */
    public String getReferenceGenome() {
        return referenceGenome;
    }

    /**
     * @param ReferenceGenome the ReferenceGenome to set
     */
    public void setReferenceGenome(String ReferenceGenome) {
        this.referenceGenome = ReferenceGenome;
    }

    /**
     * @return the baselinePercent
     */
    public int getBaselinePercent() {
        return baselinePercent;
    }

    /**
     * @param baselinePercent the baselinePercent to set
     */
    public void setBaselinePercent(int baselinePercent) {
        this.baselinePercent = baselinePercent;
    }
}
