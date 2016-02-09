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
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.ArrayUtils;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;
import static org.apache.commons.math3.stat.inference.TestUtils.pairedTTest;



/**
 *  process isomiR dispersion data output from the ParseSAMForMiRNAsStep
 * 
 *   Input is a tab delimited file 
 * 
 * @author sr
 */

public class StepAnalyzeIsomiRDispersion extends NGSStep{
    
    static Logger                       logger                      = LogManager.getLogger();
    
    private static final String         infileExtension             = ".disp.summary.tsv";
    private static final String         dispersionResultsExtension  = ".disp.test.tsv";
    
    private ArrayList<String>           sourceList;
    private ArrayList<String>           conditionList;
    private ConditionSet                pairedSet;

    private List<MiRNAFeature>          miRNAList                   = new ArrayList<>();
    private ArrayList<IsomiRSet>        isomiRList;
    private ArrayList<MiRNADispAnalysisResult> miRNAdispAnalysisResList;
    
    
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepAnalyzeIsomiRDispersion(StepInputData sid){
        stepInputData = sid;
    }
    
    @Override
    public void execute(){
        this.setPaths();
        /*
            isomiRDispersionAnalysisParams.put("pvalue", this.getAnalyzeIsomiRDispPVal());
            
        */
        try{
            stepInputData.verifyInputData();            
        }
        catch(IOException exIO){
            logger.info("exception parsing InputData" + exIO);
        }

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
        
        try{


            String isomiRDispFile = inFolder + FILESEPARATOR + stepInputData.getProjectID() + infileExtension;
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
        }
        
        logger.info("read " + isomiRList.size() + " entries");
        
        
        try{
            this.loadMiRBaseData((String) stepInputData.getStepParams().get("host"), (String) stepInputData.getStepParams().get("miRBaseHostGFFFile"));
        }
        catch(IOException ex){
            logger.error("error reading miRBase reference file <" + (String) stepInputData.getStepParams().get("miRBaseHostGFFFile") + ">\n" + ex.toString());
        }
        
        
        logger.info("testing...\n");
        
        for(MiRNAFeature miRFeature: miRNAList){
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
        
        
        String dispersionFile   = outFolder + FILESEPARATOR + stepInputData.getProjectID() + dispersionResultsExtension;
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
        String matureFAFile = new File(miRBaseGFFFile).getParent() + FILESEPARATOR + "mature.fa";
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
                String seq = miRBaseSeq.get(id);
                if(seq != null) 
                    this.miRNAList.add(new MiRNAFeature(name, chr, startPos, endPos, strand, id, parent, seq));
                else
                    logger.warn("no sequence found for entry <" + id + ">. Skipping");
            }
        brMiR.close();
        logger.info("read " + miRNAList.size() + "miRNA entries");
        
    }
    
    
    
    @Override
    public void verifyOutputData(){
        
    }
}
