/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;


/**
 *  Adapter Trimming Step for single end reads
 *  Performs adapter trimming via a call to Trimmomatic.
 * 
 *   Input is a raw FASTQ file
 *   Output is a trimmed FASTQ file
 * 
 * 
 * @author sr
 */

public class StepSingleTrimAdapters extends NGSStep{
    
    static Logger                       logger = LogManager.getLogger();
    
    public static final String          STEP_ID_STRING          = "ADAPTER_TRIMMING_SINGLE:";
    private static final String         ID_SOFTWARE             = "adapter_software:";
    private static final String         ID_ADAPTOR_FILE         = "adapter_file:";
    private static final String         ID_MISMATCHES           = "no_of_mismatches:";
    private static final String         ID_MIN_ALIGN_SCORE      = "min_align_score:";
    private static final String         ID_MIN_AVGREAD_QUAL     = "min_avg_read_qual:";
    private static final String         ID_THREADS              = "no_of_threads:";

    private static final String         INFILE_EXTENSION        = ".fastq";
    private static final String         OUTFILE_EXTENSION       = ".trim.fastq";
    
    private String                      trimSoftware            = "";
    private String                      adapterFile             = "";
    private int                         noOfThreads             = 4;
    private int                         noOfMismatches          = 2;
    private int                         minAlignScore           = 7;
    private int                         minAvgReadQuality       = 30;
    
    
    /**
     * @return the STEP_ID_STRING
     */
    public static String getSTEP_ID_STRING() {
        return STEP_ID_STRING;
    }
    
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepSingleTrimAdapters(StepInputData sid){
        try{
            stepInputData = sid;
            stepInputData.verifyInputData();
            
        }
        catch(IOException exIO){
            
        }
    }
    
    
    
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        
        if(configData.get(ID_SOFTWARE)==null) {
            throw new NullPointerException("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_ADAPTOR_FILE)==null) {
            throw new NullPointerException("<" + configData.get(ID_ADAPTOR_FILE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MISMATCHES)==null) {
            throw new NullPointerException("<" + configData.get(ID_MISMATCHES) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIN_ALIGN_SCORE)==null) {
            throw new NullPointerException("<" + configData.get(ID_MIN_ALIGN_SCORE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIN_AVGREAD_QUAL)==null) {
            throw new NullPointerException("<" + configData.get(ID_MIN_AVGREAD_QUAL) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_THREADS)==null) {
            throw new NullPointerException("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
        }
        
        
        try{
            Integer.parseInt((String) configData.get(ID_MISMATCHES));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_MISMATCHES + " <" + configData.get(ID_MISMATCHES) + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_MISMATCHES)) <= 0){
            throw new IllegalArgumentException(ID_MISMATCHES + " <" + (String) configData.get(ID_MISMATCHES) + "> must be positive");
        }
        this.setNoOfMismatches(Integer.parseInt((String) configData.get(ID_MISMATCHES)));
        
        
        try{
            Integer.parseInt((String) configData.get(ID_MIN_ALIGN_SCORE));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_MIN_ALIGN_SCORE + " <" + configData.get(ID_MIN_ALIGN_SCORE) + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_MIN_ALIGN_SCORE)) <= 0){
            throw new IllegalArgumentException(ID_MIN_ALIGN_SCORE + " <" + (String) configData.get(ID_MIN_ALIGN_SCORE) + "> must be positive");
        }
        this.setMinAlignScore(Integer.parseInt((String) configData.get(ID_MIN_ALIGN_SCORE)));
        
        
        try{
            Integer.parseInt((String) configData.get(ID_MIN_AVGREAD_QUAL));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_MIN_AVGREAD_QUAL + " <" + configData.get(ID_MIN_AVGREAD_QUAL) + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_MIN_AVGREAD_QUAL)) <= 0){
            throw new IllegalArgumentException(ID_MIN_AVGREAD_QUAL + " <" + (String) configData.get(ID_MIN_AVGREAD_QUAL) + "> must be positive");
        }
        this.setMinAvgReadQuality(Integer.parseInt((String) configData.get(ID_MIN_AVGREAD_QUAL)));
        
        
        try{
            Integer.parseInt((String) configData.get(ID_THREADS));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_THREADS)) <= 0){
            throw new IllegalArgumentException(ID_THREADS + " <" + (String) configData.get(ID_THREADS) + "> must be positive");
        }
        this.setNoOfThreads(Integer.parseInt((String) configData.get(ID_THREADS)));

        logger.info("passed");
    }
    
    
    
    
    
    
    @Override
    public void execute(){
        
        this.setPaths();
        
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                ArrayList<String> cmd = new ArrayList<>();
                cmd.add("java -jar");
                cmd.add(this.getTrimSoftware());
                cmd.add("SE");
                cmd.add("-phred64");
                cmd.add("-threads " + this.getNoOfThreads());
                cmd.add(inFolder + FILESEPARATOR + sampleData.getFastqFile1());

                Boolean f = new File(outFolder).mkdir();       
                if (f) logger.info("created output folder <" + outFolder + "> for results" );

                cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(INFILE_EXTENSION, OUTFILE_EXTENSION));
                cmd.add("ILLUMINACLIP:" + this.getAdapterFile()
                    + ":" + this.getNoOfMismatches()
                    + ":30"
                    + ":" + this.getMinAlignScore()
                );
                cmd.add("AVGQUAL:" + this.getMinAvgReadQuality());

                /*
                java -jar Trimmomatic-0.33/trimmomatic-0.33.jar 
                SE 
                -phred64	
                <input file fastq file>
                <trimmed output fastq file> 	
                ILLUMINACLIP:Trimmomatic-0.33/adapters/TruSeqE-SE.fa:2:30:7 2:30:7 (Mismatches:N/A:Alignment Score)
                AVGQUAL (drop read if average quality is below this threshold)
                -trimlog (write log)
                -threads
                */      

                String cmdTrimAdapters = StringUtils.join(cmd, " ");
                cmdTrimAdapters = cmdTrimAdapters.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR);
                logger.info("Adapter Trim command:\t" + cmdTrimAdapters);

                Runtime rt = Runtime.getRuntime();
                Process proc = rt.exec(cmdTrimAdapters);
                BufferedReader brStdin  = new BufferedReader(new InputStreamReader(proc.getInputStream()));
                BufferedReader brStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
                
                    String line = null;
                    logger.info("<OUTPUT>");
                    while ( (line = brStdin.readLine()) != null)
                        logger.info(line);
                    logger.info("</OUTPUT>");

                    logger.info("<ERROR>");
                    while ( (line = brStdErr.readLine()) != null)
                        logger.info(line);
                    logger.info("</ERROR>");
                
                
                int exitVal = proc.waitFor();            
                logger.info("Process exitValue: " + exitVal);   
                
                brStdin.close();
                brStdErr.close();
            }
            catch(IOException|InterruptedException ex){
                logger.error("error executing AdapterTrimming command\n" + ex.toString());
            }
        }
        
        
    }
    
    
    
            
    @Override
    public void verifyInputData() throws IOException{
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            String fastqFile1 = (String)sampleData.getFastqFile1();
            String fastqFile2 = (String)sampleData.getFastqFile2();
            
            if (fastqFile1==null) throw new IOException("no Fastq1 file specified");
            
            if ((new File(fastqFile1)).exists()==false){
                throw new IOException("AdapterTrimming: fastq File1 <" 
                  + fastqFile1 + "> does not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false){
                throw new IOException("AdapterTrimming: incorrect file extension for input file <" 
                  + fastqFile1 + ">. " 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
                        
        }
    }


    
    @Override
    public void verifyOutputData(){
        
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
        HashMap paramData = new HashMap();        
        
        paramData.put(StepSingleTrimAdapters.ID_SOFTWARE, "trimmomatic/Trimmomatic-0.33/trimmomatic-0.33.jar");
        paramData.put(StepSingleTrimAdapters.ID_ADAPTOR_FILE, "/Users/simonray/NetBeansProjects/SmallRNAPipeline/test/TruSeqE-SE.fa");
        paramData.put(StepSingleTrimAdapters.ID_MISMATCHES, 2);
        paramData.put(StepSingleTrimAdapters.ID_MIN_ALIGN_SCORE, 7);
        paramData.put(StepSingleTrimAdapters.ID_MIN_AVGREAD_QUAL, 30);        
        paramData.put(StepSingleTrimAdapters.ID_THREADS, 4);
        configData.put(STEP_ID_STRING, paramData);

        
        return configData;
    }
    
    
    

    /**
     * @return the trimSoftware
     */
    public String getTrimSoftware() {
        return trimSoftware;
    }

    /**
     * @param trimSoftware the trimSoftware to set
     */
    public void setTrimSoftware(String trimSoftware) {
        this.trimSoftware = trimSoftware;
    }

    /**
     * @return the noOfThreads
     */
    public int getNoOfThreads() {
        return noOfThreads;
    }

    /**
     * @param noOfThreads the noOfThreads to set
     */
    public void setNoOfThreads(int noOfThreads) {
        this.noOfThreads = noOfThreads;
    }

    /**
     * @return the adapterFile
     */
    public String getAdapterFile() {
        return adapterFile;
    }

    /**
     * @param adapterFile the adapterFile to set
     */
    public void setAdapterFile(String adapterFile) {
        this.adapterFile = adapterFile;
    }

    /**
     * @return the noOfMismatches
     */
    public int getNoOfMismatches() {
        return noOfMismatches;
    }

    /**
     * @param noOfMismatches the noOfMismatches to set
     */
    public void setNoOfMismatches(int noOfMismatches) {
        this.noOfMismatches = noOfMismatches;
    }

    /**
     * @return the minAlignScore
     */
    public int getMinAlignScore() {
        return minAlignScore;
    }

    /**
     * @param minAlignScore the minAlignScore to set
     */
    public void setMinAlignScore(int minAlignScore) {
        this.minAlignScore = minAlignScore;
    }

    /**
     * @return the minAvgReadQuality
     */
    public int getMinAvgReadQuality() {
        return minAvgReadQuality;
    }

    /**
     * @param minAvgReadQuality the minAvgReadQuality to set
     */
    public void setMinAvgReadQuality(int minAvgReadQuality) {
        this.minAvgReadQuality = minAvgReadQuality;
    }

}
