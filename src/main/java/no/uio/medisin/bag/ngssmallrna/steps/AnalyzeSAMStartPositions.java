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
import no.uio.medisin.bag.ngssmallrna.pipeline.GenomeFeature;
import no.uio.medisin.bag.ngssmallrna.pipeline.GenomeFeatureSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.IsomiRSet;
import no.uio.medisin.bag.ngssmallrna.pipeline.MiRNAFeature;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 *  parse SAM file to extract the distribution of start positions.
 *  This is to try and determine whether reads correspond to genuine small RNAs
 *  or are a consequence of mRNA degradation. If the latter, then we would find 
 *  most reads are located within coding regions (I think)
 * 
 *   Input is a SAM file
 * 
 * @author sr
 */

public class AnalyzeSAMStartPositions extends NGSStep{
    
    static Logger                       logger                      = LogManager.getLogger();
    static  String                      FileSeparator               = System.getProperty("file.separator");
    
    private static final String         inFolder                    = "bowtie_genome_mapped";
    private static final String         startAnalysisOutFolder      = "start_pos_analysis";
    
    
    private static final String         infileExtension             = ".trim.clp.gen.sam";
    private static final String         readsInExonExtension      = ".trim.clp.gen.iso_summary.tsv";
    private static final String         readsInNonCodingExtension       = ".trim.clp.gen.iso_pretty.tsv";
    
    
    
    private StepInputData               stepInputData;
    private StepResultData              stepResultData;
    
    private GenomeFeatureSet            genomeFeatureSet;
    private List<HashMap>               codingHits;
    private List<HashMap>               nonCodingHits;
    private ArrayList<IsomiRSet>        isomiRList;
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public AnalyzeSAMStartPositions(StepInputData sid){
        genomeFeatureSet = new GenomeFeatureSet();
        stepInputData = sid;
    }
    
    @Override
    public void execute(){
        /*
            analyzeSAMStartPositionsParams.put("bleed", this.getSamParseStartPosBleed());
            analyzeSAMStartPositionsParams.put("feature_types", this.getSamParseFeatureTypes());
            analyzeSAMStartPositionsParams.put("host", this.getBowtieMappingReferenceGenome());
            analyzeSAMStartPositionsParams.put("genomeReferenceGFFFile", this.getGenomeAnnotationGFF());
        
        */
        try{
            stepInputData.verifyInputData();            
        }
        catch(IOException exIO){
            logger.info("exception parsing InputData" + exIO);
        }
    
        
        try{
            this.genomeFeatureSet.loadGenomeGFFData((String) stepInputData.getStepParams().get("genomeReferenceGFFFile"));
        }
        catch(IOException ex){
            logger.error("error reading genome reference GFF file <" + (String) stepInputData.getStepParams().get("genomeReferenceGFFFile") + ">\n" + ex.toString());
        }
        ArrayList<String> x = ((ArrayList<String>)stepInputData.getStepParams().get("feature_types"));
        String[] featureTypes = x.toArray(new String[x.size()]);

        String samInputFile = null;
        String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
        String startAnalysisOutputFolder = pathToData + FileSeparator + startAnalysisOutFolder;
        startAnalysisOutputFolder = startAnalysisOutputFolder.replace(FileSeparator + FileSeparator, FileSeparator).trim();
        Boolean fA = new File(startAnalysisOutputFolder).mkdir();       
        if (fA) logger.info("created output folder <" + startAnalysisOutputFolder + "> for results" );
        
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        codingHits      = new ArrayList<>();
        nonCodingHits   = new ArrayList<>();
        while (itSD.hasNext()){
            try{
                
                
                int bleed = (int) stepInputData.getStepParams().get("bleed");
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                
                samInputFile = pathToData + FileSeparator + inFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", infileExtension);
                logger.info(sampleData.getDataFile().replace(".fastq", infileExtension));
                int matchCount5 = 0;
                int matchCount3 = 0;
                int codingCounts = 0;
                int nonCodingCounts = 0;
                int preMatchCount5 = 0;
                int preMatchCount3 = 0;
                int totalCounts = 0;
                int lineCount = 0;
                String samLine = null;
                BufferedReader brSAM = new BufferedReader(new FileReader(new File(samInputFile)));
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
                            String name = samLine.split("\t")[0];
                            String chr = samLine.split("\t")[2].trim();
                            String mdString = samLine.split("\t")[12];
                            
                            GenomeFeature featureHit = this.genomeFeatureSet.doesReadOverlapFeature(startPos, endPos, chr, strand, bleed, featureTypes);
                            if (featureHit == null){
                                
                                HashMap ncRead = new HashMap();
                                ncRead.put("start", startPos);
                                ncRead.put("end", endPos);
                                ncRead.put("name", name);
                                ncRead.put("chr", chr);
                                ncRead.put("strand", strand);
                                ncRead.put("counts", Integer.parseInt(samLine.split("\t")[0].split("-")[1]));
                                nonCodingHits.add(ncRead);
                                nonCodingCounts += Integer.parseInt(samLine.split("\t")[0].split("-")[1]);
                                
                            }
                            else{
                                
                                HashMap cRead = new HashMap();
                                cRead.put("start", featureHit.getStartPos() - startPos);
                                cRead.put("end", featureHit.getStartPos() - endPos);
                                cRead.put("name", name);
                                cRead.put("chr", chr);
                                cRead.put("strand", strand);
                                cRead.put("counts", Integer.parseInt(samLine.split("\t")[0].split("-")[1]));
                                codingHits.add(cRead);
                                codingCounts += Integer.parseInt(samLine.split("\t")[0].split("-")[1]);
                                
                            }
                            
                            
                            
                            
                        }
                        lineCount++;
                        if((lineCount % 10000) == 0)
                            logger.info(lineCount);
                        
                    }
                    logger.info("  total reads              = " + totalCounts);
                    logger.info("  total coding reads       = " + codingCounts);
                    logger.info("  total non coding reads   = " + nonCodingCounts);
                    Double minCounts = (double) totalCounts /100000.0;
                    logger.info((matchCount5 + matchCount3) + " reads (" + matchCount5 + " 5'" + "/" + matchCount3 + " 3' ) were mapped");
                    String  codingReadsDetailsFile = startAnalysisOutputFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", readsInExonExtension);
                    String  nonCodingReadsDetailsFile  = startAnalysisOutputFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", readsInNonCodingExtension);
                    
            
                    
                    try{
                    logger.info("  write coding information");
                        BufferedWriter brCoding = new BufferedWriter(new FileWriter(new File(codingReadsDetailsFile)));
                            brCoding.write("name\tstart\tend\tchr\tstrand\tcounts\n");
                            for(HashMap codingHit: codingHits){
                                brCoding.write(codingHit.get("name") + "\t" + codingHit.get("start") + "\t" + codingHit.get("end") 
                                  + "\t" + codingHit.get("chr") + "\t" + codingHit.get("strand") + "\t" + codingHit.get("counts") + "\n");                        
                            }
                        brCoding.close();
                    }
                    catch(IOException exIO){
                        logger.info("error writing coding File <" + codingReadsDetailsFile + ">\n" + exIO);
                    }

                    try{
                        logger.info("  write non-coding information");                    
                        BufferedWriter brNonCoding  = new BufferedWriter(new FileWriter(new File(nonCodingReadsDetailsFile)));
                            brNonCoding.write("name\tstart\tend\tchr\tstrand\tcounts\n");
                            for(HashMap nonCodingHit: nonCodingHits){
                                brNonCoding.write(nonCodingHit.get("name") + "\t" + nonCodingHit.get("start") + "\t" + nonCodingHit.get("end") 
                                  + "\t" + nonCodingHit.get("chr") + "\t" + nonCodingHit.get("strand") + "\t" + nonCodingHit.get("counts") + "\n");                        
                            }
                        brNonCoding.close();
                    }
                    catch(IOException exIO){
                        logger.info("error writing non-coding File <" + nonCodingReadsDetailsFile + ">\n" + exIO);
                    }
                    
                brSAM.close();
                logger.info("  done\n");
            }
            catch(IOException exIO){
                logger.info("error reading non-coding File <" + samInputFile + ">\n" + exIO);
            }
                
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
    
    
    
    
    @Override
    public void outputResultData(){
        
    }
}
