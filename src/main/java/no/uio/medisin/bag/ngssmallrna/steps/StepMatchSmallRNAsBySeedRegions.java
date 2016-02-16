/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import com.abahgat.suffixtree.GeneralizedSuffixTree;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.MirFeatureSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;




/**
 *  Adapter Trimming Step
 *  Unzip FASTQ files using pigz.
 * 
 *   Input is a zipped FASTQ file
 *   Output is a unzipped FASTQ file
 * 
 * 
 * @author sr
 */

public class StepMatchSmallRNAsBySeedRegions extends NGSStep implements NGSBase{
    
    static Logger                       logger = LogManager.getLogger();
    
    public static final String          STEP_ID_STRING          = "MatchSmallRNAsBySeedRegions";
    private static final String         ID_MIRBASE_VERSION      = "mirbaseVersion";
    private static final String         ID_REF_GENOME           = "host";
    private static final String         ID_THREADS              = "noOfThreads";
    
    private static final String         INFILE_EXTENSION        = ".fastq.gz";
    private static final String         OUTFILE_EXTENSION       = ".fastq";
    
    MirFeatureSet                       mirBaseSet              = new MirFeatureSet();  
    
    private int                         noOfThreads             = 4;
    private String                      unzipSoftware           = "";
    private int                         miRBaseRelease          = 20;
    private String                      referenceGenome         = "";

    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepMatchSmallRNAsBySeedRegions(StepInputData sid){
       stepInputData = sid;
    }
    
    
    
    /**
     * This parses out the hashmap containing the run parameters for this step
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        if(configData.get(ID_THREADS)==null) {
            logger.error("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
        }
        
        
        try{
            this.setNoOfThreads((Integer)configData.get(ID_THREADS));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_THREADS + " <" + configData.get(ID_THREADS) + "> is not an integer");
            throw new NumberFormatException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> is not an integer");
        }
        
        if (this.getNoOfThreads() <= 0){
            logger.error(ID_THREADS + " <" + configData.get(ID_THREADS) + "> must be positive");
            throw new IllegalArgumentException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> must be positive");
        }
        
 
        logger.info("passed");
    }
    
    
    
    /**
     * unzip specified file types
     * 
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException{
        /*
            sample unzip command       
            pigz -p 4 -d /data/ngsdata/project1/sra_data.fastq.gz     
        */
        logger.info(STEP_ID_STRING + ": execute step");        
        
        GeneralizedSuffixTree in = new GeneralizedSuffixTree();

        String word = "cacao";
        in.put(word, 0);

        String gffFileMirBase = this.cleanPath(stepInputData.getDataLocations().getMirbaseFolder() 
                + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3");
        String faFileMirBase = gffFileMirBase.replace("gff3", "mature.fa");
        mirBaseSet.loadMiRBaseData(this.getReferenceGenome(), gffFileMirBase, faFileMirBase);
        
        // cycle through each entry and search for conserved family match
        
        logger.info(STEP_ID_STRING + ": completed");
    }
    
    
    
            
    /**
     * this should be called prior to executing the step.
     * check unzip software exists and input files are available
     * 
     * @throws IOException
     */
    @Override
    public void verifyInputData() throws IOException{
        
        logger.info("verify input data");        
        this.setPaths();
        
        if(new File(this.getUnzipSoftware()).exists() == false){
            logger.error("unzip software not found at location < " + this.getUnzipSoftware() +">");
            throw new IOException("unzip software not found at location < " + this.getUnzipSoftware() +">");
        }
                
                
                            
        // check the data files
        String fastqFile1in = "";
        String fastqFile1out = "";
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            
            //Fastq 1
            if (sampleData.getFastqFile1()==null) {
                logger.error("no Fastq1 file specified");
                throw new IOException("no Fastq1 file specified");
            }
            
            fastqFile1out = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1());
            fastqFile1in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION));
            String fastqFile1 = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION);
            if (new File(fastqFile1in).exists()==false && new File(fastqFile1out).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq 1 Files <" + fastqFile1in + "> & <" + fastqFile1out + "> do not exist");
                throw new IOException(STEP_ID_STRING + ": fastq 2 Files <" + fastqFile1in + "> & <" + fastqFile1out + "> do not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.info(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1 + ">. should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
            
            //Fastq 2
            if (sampleData.getFastqFile2()==null) continue;
            String fastqFile2in = "";
            String fastqFile2out = "";
            fastqFile2in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(".fastq", INFILE_EXTENSION));
            fastqFile2out = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile2());

            
            if ((new File(fastqFile2in)).exists()==false && (new File(fastqFile2in)).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq 2 Files <" + fastqFile2in + "> & <" + fastqFile2out + "> do not exist");
                throw new IOException(STEP_ID_STRING + ": fastq 2 Files <" + fastqFile2in + "> & <" + fastqFile2out + "> do not exist");
            }
            if (fastqFile2in.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                  + fastqFile2in + ">. should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                  + fastqFile2in + ">. \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
                        
            
        }

    }
    
    
    
    @Override
    public void verifyOutputData(){
        logger.info("no output verification required");
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
        
        HashMap<String, Object> configData = new HashMap();
        
        configData.put(ID_REF_GENOME, "hsa");
        configData.put(ID_MIRBASE_VERSION, 20);
        configData.put(ID_THREADS, 4);
        
        return configData;
        
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
     * @return the unzipSoftware
     */
    public String getUnzipSoftware() {
        return unzipSoftware;
    }

    /**
     * @param unzipSoftware the unzipSoftware to set
     */
    public void setUnzipSoftware(String unzipSoftware) {
        this.unzipSoftware = unzipSoftware;
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
     * @return the referenceGenome
     */
    public String getReferenceGenome() {
        return referenceGenome;
    }

    /**
     * @param referenceGenome the referenceGenome to set
     */
    public void setReferenceGenome(String referenceGenome) {
        this.referenceGenome = referenceGenome;
    }
    
    
    
    
    
}
