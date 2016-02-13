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
import no.uio.medisin.bag.ngssmallrna.pipeline.ConditionEntry;
import no.uio.medisin.bag.ngssmallrna.pipeline.ConditionSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.MiRNADispAnalysisResult;
import no.uio.medisin.bag.ngssmallrna.pipeline.IsomiRSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.MiRNAFeature;
import no.uio.medisin.bag.ngssmallrna.pipeline.MirFeatureSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.ArrayUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;



/**
 *  process isomiR dispersion data output from the ParseSAMForMiRNAsStep
 * 
 *   Input is a tab delimited file 
 * 
 * @author sr
 */

public class StepAnalyzeIsomiRDispersion extends NGSStep implements NGSBase{
    
    static Logger                   logger                          = LogManager.getLogger();
    public  static final String     STEP_ID_STRING                  = "AnalyzeIsomiRDispersion";
    private static final String     ID_REF_GENOME                   = "host";
    private static final String     ID_MIRBASE_VERSION              = "mirbaseVersion";
    private static final String     ID_PVALUE                       = "pValue:";
    
    private static final String     INFILE_EXTENSION                = ".disp.summary.tsv";
    private static final String     DISPERSION_RESULTS_EXTENSION    = ".disp.test.tsv";
    
    private ArrayList<String>       sourceList;
    private ArrayList<String>       conditionList;
    private ConditionSet            pairedSet;

    //private List<MiRNAFeature>      miRNAList                       = new ArrayList<>();
    MirFeatureSet                   mirBaseSet                      = new MirFeatureSet();           
    private ArrayList<IsomiRSet>    isomiRList;
    private ArrayList<MiRNADispAnalysisResult> miRNAdispAnalysisResList;
    
    private double                  pValue                          = 0.05;
    private String                  referenceGenome                 = "";
    private int                     miRBaseRelease                  = 20;
    
    
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepAnalyzeIsomiRDispersion(StepInputData sid){
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
        
        if(configData.get(ID_PVALUE)==null) {
            throw new NullPointerException("<" + ID_PVALUE + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIRBASE_VERSION)==null) {
            throw new NullPointerException("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
        }
        

      
        try{
            this.setMiRBaseRelease((Integer) configData.get(ID_MIRBASE_VERSION));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_MIRBASE_VERSION + " <" + configData.get(ID_MIRBASE_VERSION) + "> is not an integer");
        }        
        if (Integer.parseInt((String) configData.get(ID_MIRBASE_VERSION)) <= 0){
            throw new IllegalArgumentException(ID_MIRBASE_VERSION + " <" + configData.get(ID_MIRBASE_VERSION) + "> must be positive integer");
        }
        

        
        try{
            this.setpValue((Double) configData.get(ID_PVALUE));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_PVALUE + " <" + configData.get(ID_PVALUE) + "> is not an integer");
        }        
        if (Double.parseDouble((String) configData.get(ID_PVALUE)) <= 0 || Double.parseDouble((String) configData.get(ID_PVALUE)) > 1.0){
            throw new IllegalArgumentException(ID_PVALUE + " <" + configData.get(ID_PVALUE) + "> must be an float between 0.0 and 1.0");
        }
        

        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }
        

        logger.info("passed");
    }
    
    
    
    @Override
    public void execute() throws IOException{
        this.setPaths();

        stepInputData.verifyInputData();            

        sourceList      = new ArrayList<>();
        conditionList   = new ArrayList<>();
        pairedSet       = new ConditionSet();

        
        /*
            generate paired data points
        */
        for(SampleDataEntry dataEntry: stepInputData.getSampleData()){
            
            if(sourceList.contains(dataEntry.getDataSource())==false){
                pairedSet.addEntry(dataEntry.getDataSource(), dataEntry.getFastqFile1().replace(".fastq", ""));
                sourceList.add(dataEntry.getDataSource());
            }
            else{
                pairedSet.updateEntry(dataEntry.getDataSource(), dataEntry.getFastqFile1().replace(".fastq", ""));
            }
            
            if(conditionList.contains(dataEntry.getCondition())==false){
                conditionList.add(dataEntry.getCondition());
            }
            
        }
        
        logger.info("data set contains " + sourceList.size() + " unique sources and " + conditionList.size() + " unique conditions");
        int i = pairedSet.removeUnpairedEntries();
        if(i == 1)
            logger.info("removed " + i +  " paired entry");
        else{
            if(i > 1)
                logger.info("removed " + i +  " paired entries");            
        }
        Boolean fA = new File(outFolder).mkdir();       
        if (fA) logger.info("created output folder <" + outFolder + "> for results" );
        isomiRList = new ArrayList<>();
        miRNAdispAnalysisResList = new ArrayList<>();
        String isomiRDispFile = "";
        try{


            isomiRDispFile = inFolder + FILESEPARATOR + stepInputData.getProjectID() + INFILE_EXTENSION;
            logger.info("reading " + isomiRDispFile);

            String isoDispLine = null;
            BufferedReader brIDF = new BufferedReader(new FileReader(new File(isomiRDispFile)));
                String headerLine = brIDF.readLine();
                while((isoDispLine=brIDF.readLine())!= null){

                    isomiRList.add(new IsomiRSet(isoDispLine));

                }

            logger.info("  done");


        }
        catch(IOException ex){
            logger.error("error reading isomiR dispersion file\n" + ex.toString());
            throw new IOException(STEP_ID_STRING + ": error reading isomiR dispersion file <" + isomiRDispFile + ">");
        }
        
        logger.info("read " + isomiRList.size() + " entries");
        
        String gffFileMirBase = stepInputData.getDataLocations().getMirbaseFolder() + FILESEPARATOR + this.getMiRBaseRelease() + this.getReferenceGenome() + ".gff3";
        String faFileMirBase = gffFileMirBase.replace("gff3", "fasta");
        mirBaseSet.loadMiRBaseData(gffFileMirBase, faFileMirBase, this.getReferenceGenome());
        
        
        logger.info("testing...\n");
        
        for(MiRNAFeature miRFeature: this.mirBaseSet.getMiRBaseMiRNAList()){
            logger.info(miRFeature.getName());
            ArrayList<IsomiRSet>  localIsomiRList = new ArrayList();
            for(IsomiRSet isomiREntry: isomiRList){
                if(isomiREntry.getMimatID().equals(miRFeature.getMimatID())){
                    localIsomiRList.add(isomiREntry); 
                }
            }
            
            
            if(localIsomiRList.size() > 3){
                MiRNADispAnalysisResult miRNADispAnalRes = this.analyzeIsomiRSet(miRFeature.getMimatID(), localIsomiRList);

                miRNAdispAnalysisResList.add(miRNADispAnalRes);
            }
            
            
        }
        
        
        String dispersionFile   = outFolder + FILESEPARATOR + stepInputData.getProjectID() + DISPERSION_RESULTS_EXTENSION;
        logger.info("write dispersions to file <" + dispersionFile + ">");
        try{
            BufferedWriter bwDp = new BufferedWriter(new FileWriter(new File(dispersionFile)));
                bwDp.write(MiRNADispAnalysisResult.printHeaderLine());
                for(MiRNADispAnalysisResult miRNAdispAnalysisRes: miRNAdispAnalysisResList){
                    bwDp.write(miRNAdispAnalysisRes.printSummary());                
                }
            bwDp.close();
        }
        catch(IOException exIO){
            logger.info("error writing isomiR dispersion File <" + dispersionFile + ">\n" + exIO);
            throw new IOException(STEP_ID_STRING + ": error writing isomiR dispersion File <" + dispersionFile + ">");
        }
                    
        
    }
    
    
    
    
    /**
     * perform some kind of analysis on this isomiR set
     * This version assumes the conditions exist in pairs, i.e. there are
     * two NGS files for each condition (patient). Tumour/non-tumour
     * this will obviously require further development but is merely
     * used here for a test of the basic idea
     * 
     * @param  mimatID
     * @param  localIsomiRList ArrayList(IsomiRSet)
     * @return Point2D.Double
     */
    public MiRNADispAnalysisResult analyzeIsomiRSet(String mimatID, ArrayList<IsomiRSet>  localIsomiRList){
        List<Double> sample1 = new ArrayList();
        List<Double> sample2 = new ArrayList();
        String note = "";

        logger.info(pairedSet.printPairedList());
        for(ConditionEntry conditionEntry: pairedSet.getPairedList()){
            
            // i want to get paired entries, one for each source 
            double s1 = 0.0;
            double s2 = 0.0;
            boolean b1 = false;
            boolean b2 = false;
            String note1= "";
            for(IsomiRSet isomiRset: localIsomiRList){
                
                if(isomiRset.getRunID().equals(conditionEntry.getDataset1())){
                    s1 = isomiRset.getTotal3pDist();
                    b1 = true;
                }
                else{
                    if(isomiRset.getRunID().equals(conditionEntry.getDataset2())){
                        s2 = isomiRset.getTotal3pDist();
                        b2 = true;
                        if (b1) {
                            note1 = isomiRset.getNote();
                        }
                    }
                }
            }
            if(b1 && b2){
                if(s1 != 0.0 && s2 != 0.0)
                {
                    sample1.add(s1);
                    sample2.add(s2);            
                    note = note.concat(note1 + ":");
                }
            }
            
            
        }

        for(double s: sample1){
            logger.info(s);
        }
        for(double s: sample2){
            logger.info(s);
        }
        double pValue = 1.0;
        if((Integer)sample1.size()>3){
            
            pValue =  pairedTTest(
              ArrayUtils.toPrimitive(sample1.toArray(new Double[sample1.size()])), 
              ArrayUtils.toPrimitive(sample2.toArray(new Double[sample2.size()])));

        }
        Double span = ((Integer)sample1.size()).doubleValue() / ((Integer) pairedSet.getNumberOfEntries()).doubleValue();
        logger.info("  " + pValue + "\t\t" + span);
        return new MiRNADispAnalysisResult(mimatID, span, pValue, note);
       
    }
    
    
    /**
     * Verify Input Data 
     * 
     * @throws IOException
     */        
    @Override
    public void verifyInputData() throws IOException{
        
        logger.info("verify input data");        
        this.setPaths();
                
        String isomiRDispFile = inFolder + FILESEPARATOR + stepInputData.getProjectID() + INFILE_EXTENSION;
        String gffFileMirBase = stepInputData.getDataLocations().getMirbaseFolder() + FILESEPARATOR + this.getMiRBaseRelease() + this.getReferenceGenome() + ".gff3";
        String faFileMirBase = gffFileMirBase.replace("gff3", "fasta");

        if(new File(isomiRDispFile).exists() == false){
            throw new IOException(STEP_ID_STRING + ": isomiR dispersion not found at location < " + isomiRDispFile +">");
        }
        
        if(new File(gffFileMirBase).exists() == false){
            throw new IOException(STEP_ID_STRING + ": miRBase GFF file not found at location < " + gffFileMirBase +">");
        }
        
        if(new File(faFileMirBase).exists() == false){
            throw new IOException(STEP_ID_STRING + ": miRBase Fasta file not found at location < " + faFileMirBase +">");
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
        configData.put(ID_PVALUE, 2);
        configData.put(ID_MIRBASE_VERSION, 20);


        return configData;
    }


    
    
    
    @Override
    public void verifyOutputData(){
        
    }

    /**
     * @return the pValue
     */
    public double getpValue() {
        return pValue;
    }

    /**
     * @param pValue the pValue to set
     */
    public void setpValue(double pValue) {
        this.pValue = pValue;
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
}
