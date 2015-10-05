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

public class DEwithEdgeRStep extends NGSStep{
    
    static Logger                       logger                      = LogManager.getLogger();
    static  String                      FileSeparator               = System.getProperty("file.separator");
    
    private static final String         inFolder                    = "mirna_isomir_analysis";
    private static final String         deAnalysisOutFolder         = "de_analysis";
    
    
    private static final String         infileExtension             = ".trim.clp.gen.sam";
    private static final String         isomirSummaryExtension      = ".trim.clp.gen.iso_summary.tsv";
    private static final String         isomirPrettyExtension       = ".trim.clp.gen.iso_pretty.tsv";
    private static final String         miRCountsExtension          = ".trim.clp.gen.mircounts.tsv";
    
    
    
    private StepInputData               stepInputData;
    private StepResultData              stepResultData;
    

    private List<MiRNAFeature>          miRNAList                   = new ArrayList<>();
    private List<MiRNAFeature>          miRNAHitList;
    private ArrayList<IsomiRSet>        isomiRList;
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public DEwithEdgeRStep(StepInputData sid){
        stepInputData = sid;
    }
    
    @Override
    public void execute(){
        /*
            diffExpressionAnalysisParams.put("pvalue", this.getDiffExpressionPVal());
        */
        try{
            stepInputData.verifyInputData();            
        }
        catch(IOException exIO){
            logger.info("exception parsing InputData" + exIO);
        }
    
        /*
            1. read in all sample count files and merge
            2. output merged count file
            3. generate R script to perform DE using EdgeR
            4. process EdgeR output file 
        try{
            this.loadMiRBaseData((String) stepInputData.getStepParams().get("host"), (String) stepInputData.getStepParams().get("miRBaseHostGFFFile"));
        }
        catch(IOException ex){
            logger.error("error reading miRBase reference file <" + (String) stepInputData.getStepParams().get("miRBaseHostGFFFile") + ">\n" + ex.toString());
        }
        */
        
        
        String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
        String miRNAInputFolder = pathToData + FileSeparator + inFolder;
        String deAnalysisOutputFolder = pathToData + FileSeparator + deAnalysisOutFolder;
        
        deAnalysisOutputFolder = deAnalysisOutputFolder.replace(FileSeparator + FileSeparator, FileSeparator).trim();
        Boolean fA = new File(deAnalysisOutputFolder).mkdir();       
        if (fA) logger.info("created output folder <" + deAnalysisOutputFolder + "> for results" );
        
        Iterator itSD = this.stepInputData.getSampleData().iterator();

        
        logger.info("Merging Count Files");
        String[] countStrings = new String[miRNAList.size()];
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            String  miRCountsFile  = miRNAInputFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", miRCountsExtension);
            miRCountsFile = miRCountsFile.replace(FileSeparator + FileSeparator, FileSeparator).trim();
            try{
                int m=0;
                BufferedReader brmiRCounts  = new BufferedReader(new FileReader(new File(miRCountsFile)));
                    String countLine = "";
                    while((countLine = brmiRCounts.readLine()) != null){
                        countStrings[m] = countStrings[m].concat(countLine.split("\t")[1].trim());
                        m++;
                    }
                brmiRCounts.close();
            }
            catch(IOException ex){
                logger.error("error reading count files for merging <" + miRCountsFile + "> \n" + ex.toString());
            }
        }
        
        logger.info("Writing merged count files");
        String mergedCountsFile      = deAnalysisOutputFolder + FileSeparator + stepInputData.getProjectID() + ".merged.mirna_counts.tsv";    
        try{
            BufferedWriter bwMc = new BufferedWriter(new FileWriter(new File(mergedCountsFile)));
            int m=0;
            for(MiRNAFeature miR: this.miRNAList){
                bwMc.write(miR.getMimatID() + ":" + miR.getName() + "\t" + countStrings[m]);
            }
            
            bwMc.close();
        }
        catch(IOException exIO){
            logger.info("error writing merged counts File <" + mergedCountsFile + ">\n" + exIO);        
        }
        
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
        String matureFAFile = new File(miRBaseGFFFile).getParent() + FileSeparator + "mature.fa";
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
    public void outputResultData(){
        
    }
}
