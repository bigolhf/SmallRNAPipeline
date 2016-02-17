/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.GFFSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.MirFeatureSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import no.uio.medisin.bag.ngssmallrna.pipeline.TargetScanMirFamilyList;
import no.uio.medisin.bag.ngssmallrna.pipeline.TargetScanPredictedTargetList;
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
    
    private static final String         QUERYFA_EXTENSION       = ".fasta";
    private static final String         QUERYCOUNTS_EXTENSION   = ".merged.mirna_counts.tsv";
    private static final String         QUERYFEATURES_EXTENSION = ".features.tsv";
    
    MirFeatureSet                       mirBaseSet              = new MirFeatureSet(); 
    TargetScanMirFamilyList             tScanMirFamilies        = new TargetScanMirFamilyList();
    TargetScanPredictedTargetList       tscanPredictedTargets   = new TargetScanPredictedTargetList();
    
    private int                         noOfThreads             = 4;
    private String                      unzipSoftware           = "";
    private int                         miRBaseRelease          = 20;
    private String                      referenceGenome         = "";
    
    private String                      queryGenome             = "";
    private String                      queryCountDataFile      = "";
    private String                      queryFastaFile          = "";
    private String                      queryFeatureFile        = "";
    private String                      hostVersusQuerySAMFile  = "";
    
    private ArrayList<String>           querySeedSet;
    private GFFSet                      queryFeatures;
    
    
    

    

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
        
        /*
        GeneralizedSuffixTree in = new GeneralizedSuffixTree();

        String word = "cacao";
        in.put(word, 0);
        */
        
        /*
        String gffFileMirBase = this.cleanPath(stepInputData.getDataLocations().getMirbaseFolder() 
                + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3");
        String faFileMirBase = gffFileMirBase.replace("gff3", "mature.fa");
        mirBaseSet.loadMiRBaseData(this.getReferenceGenome(), gffFileMirBase, faFileMirBase);
        
        String targetScanFamily = this.cleanPath(stepInputData.getDataLocations().getTargetbaseFolder()
                + FILESEPARATOR + ReferenceDataLocations.ID_TSCAN_MIRFAMILY_FILE);
        tScanMirFamilies.loadConservedFamilyList(targetScanFamily);
        
        
        String targetScanPrediction = this.cleanPath(stepInputData.getDataLocations().getTargetbaseFolder()
                + FILESEPARATOR + ReferenceDataLocations.ID_TSCAN_MIRFAMILY_FILE);
        tscanPredictedTargets.loadPredictedTargetInfo(targetScanPrediction);
        */
        
        // need a list of tabbed delimited MSY mapped smallRNAs entries
        /*
            ultimately, this should be the following steps:
                1. load the MSY fasta file
                2. load the count data
                3. filter the fasta entries by counts        
                4. map the passed MSY sequences to the HSA genome
                5. parse the SAM file to extract the reads that have a mapped seed region
                6. 
            for now, lets start from step 5
        */
        
        this.setQueryFastaFile(this.cleanPath(this.inFolder + FILESEPARATOR + stepInputData.getProjectID() + QUERYFA_EXTENSION));
        String faLine = "";
        //querySeedSet
        /*
            for now we just take positions 2 to 8 and store these as the seed.
            in reality it is more complex because some families appear to be shifted by +1 nt 
            e.g. the seed sequence GGAGCUC is 2 to 8 in the miR-25 family, but at 3-9 in the human counterpart
        
                >grep -A 1 msy-187803 ../analyzeStartPositions500.15.30/sweden.fasta
                >msy-187803:5 187803:5 msy 187803:5        
                UUGGGGGUUGGAGAGAAGGUCCA
        
                >grep  msy-187803 ../analyzeStartPositions500.15.30/sweden.features.tsv 
                chr3	.	smallRNA	332185	332207	.	+	.	ID=msy-187803:5;Alias=msy-187803:5;Name=msy-187803:5
        
                >grep msy-187803 ../differentialExpression.500.15.30/sweden.merged.mirna_counts.tsv
                msy-187803:5|msy-187803:5	0	0	0	0	0	0	0	0	0	0    
        
        */
        queryFeatures = new GFFSet();
        this.setQueryFeatureFile(this.cleanPath(this.inFolder + FILESEPARATOR + stepInputData.getProjectID() + QUERYFEATURES_EXTENSION));
        logger.info("reading Query Feature file <" + this.getQueryFeatureFile() + ">");
        int lines = queryFeatures.readGFF(this.getQueryFeatureFile());
        logger.info("completed");
        logger.info("read " + lines + " lines");
        
        int faMatchCount = 0;
        logger.info("reading Query Fasta file <" + this.getQueryFastaFile() + ">");
        try(BufferedReader brFA = new BufferedReader(new FileReader(new File(this.getQueryFastaFile())))){
            while((faLine=brFA.readLine())!=null){
                String headerLine = faLine;
                if (headerLine.startsWith(">")==false){
                    logger.error("error reading Query Fasta file <" + this.getQueryFastaFile() + ">");
                    logger.error("error occurred on header line " + headerLine);
                    throw new IOException("error reading Query Fasta file <" + this.getQueryFastaFile() + ">");
                }
                String seqLine = brFA.readLine();
                String thisID = faLine.substring(1).split(" ")[0].trim();
                if(queryFeatures.findEntryByID(thisID)!=null){
                    queryFeatures.findEntryByID(thisID).setSequence(seqLine.trim());
                    faMatchCount++;
                }
            }
            brFA.close();
            logger.info("read " + querySeedSet.size() + "lines");
            logger.info("matched " + faMatchCount + " entries in the feature file");
        }
        catch(IOException exIO){
            logger.error("error reading Query Fasta file <" + this.getQueryFastaFile() + ">");
            logger.error("error occurred on line " + faLine);
            logger.error(exIO);
            throw new IOException("error reading Query Fasta file <" + this.getQueryFastaFile() + ">\nsee log file for details");
        }
        
        int ctMatchCount = 0;
        int dualMatchCount = 0;
        String countLine = "";
        this.setQueryCountDataFile(this.cleanPath(this.inFolder + FILESEPARATOR + stepInputData.getProjectID() + QUERYCOUNTS_EXTENSION));
        try(BufferedReader brCF = new BufferedReader(new FileReader(new File(this.getQueryCountDataFile())))){
            while((countLine=brCF.readLine())!=null){
                String thisID = faLine.split(" ")[0].trim();
                String counts[] = countLine.split("\t");
                int totalCounts = 0;
                for(int i=1; i< counts.length; i++){
                    totalCounts += Integer.parseInt(counts[i]);
                }
                int avgCounts = totalCounts / (counts.length-1);
                if(queryFeatures.findEntryByID(counts[0].trim())!=null){
                    queryFeatures.findEntryByID(thisID).addAttr("avgcounts=", Integer.toString(avgCounts));
                    ctMatchCount++;
                    if(queryFeatures.findEntryByID(counts[0].trim()).getAttr().contains("seq=")){
                        dualMatchCount++;
                    }
                }
            }
            brCF.close();
            logger.info("read " + querySeedSet.size() + "lines");
            logger.info("matched " + ctMatchCount + " entries in the feature file");
            logger.info("matched " + dualMatchCount + " entries in the feature file with 'seq=' attr");
            
        }
        catch(IOException exIO){
            logger.error("error reading Query Fasta file <" + this.getQueryFastaFile() + ">");
            logger.error("error occurred on line " + faLine);
            logger.error(exIO);
            throw new IOException("error reading Query Fasta file <" + this.getQueryFastaFile() + ">\nsee log file for details");
        }
        
        
        //
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
        
                
                
                            
        // check the data files are present
                        
            


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

    /**
     * @return the queryGenome
     */
    public String getQueryGenome() {
        return queryGenome;
    }

    /**
     * @param queryGenome the queryGenome to set
     */
    public void setQueryGenome(String queryGenome) {
        this.queryGenome = queryGenome;
    }

    /**
     * @return the queryCountDataFile
     */
    public String getQueryCountDataFile() {
        return queryCountDataFile;
    }

    /**
     * @param queryCountDataFile the queryCountDataFile to set
     */
    public void setQueryCountDataFile(String queryCountDataFile) {
        this.queryCountDataFile = queryCountDataFile;
    }

    /**
     * @return the queryFastaFile
     */
    public String getQueryFastaFile() {
        return queryFastaFile;
    }

    /**
     * @param queryFastaFile the queryFastaFile to set
     */
    public void setQueryFastaFile(String queryFastaFile) {
        this.queryFastaFile = queryFastaFile;
    }

    /**
     * @return the hostVersusQuerySAMFile
     */
    public String getHostVersusQuerySAMFile() {
        return hostVersusQuerySAMFile;
    }

    /**
     * @param hostVersusQuerySAMFile the hostVersusQuerySAMFile to set
     */
    public void setHostVersusQuerySAMFile(String hostVersusQuerySAMFile) {
        this.hostVersusQuerySAMFile = hostVersusQuerySAMFile;
    }

    /**
     * @return the queryFeatureFile
     */
    public String getQueryFeatureFile() {
        return queryFeatureFile;
    }

    /**
     * @param queryFeatureFile the queryFeatureFile to set
     */
    public void setQueryFeatureFile(String queryFeatureFile) {
        this.queryFeatureFile = queryFeatureFile;
    }
    
    
    
    
    
}
