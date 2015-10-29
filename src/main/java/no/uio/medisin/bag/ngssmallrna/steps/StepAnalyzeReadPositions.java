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
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.MappedRead;
import no.uio.medisin.bag.ngssmallrna.pipeline.SAMEntry;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 * this step performs additional analyses on a sample set of mapped reads in a SAM
 * file. It can look at features such as
 * 
 *      1. distribution of start and end positions
 *      2. distribution of mapped read length
 * 
 * 
 * @author sr
 */
public class StepAnalyzeReadPositions extends NGSStep{
    
    
    public static final int QNAME = 1;
    public static final int FLAG = 2;
    public static final int RNAME = 3;
    public static final int POS = 4;
    public static final int MAPQ = 5;
    public static final int CIGAR = 6;
    public static final int RNEXT = 7;
    public static final int PNEXT = 8;
    public static final int TLEN = 9;
    public static final int SEQ = 10;
    public static final int QUAL = 11;
    
    static Logger                       logger                      = LogManager.getLogger();
    static String                       FileSeparator               = System.getProperty("file.separator");

    
    private static final String         inFolder                    = "bowtie_genome_mapped";
    private static final String         posAnalysisFolder           = "position_analysis";

    
    private static final String         infileExtension             = ".trim.clp.gen.sam";
    private static final String         positionsExtension          = ".trim.clp.gen.pos.tsv";
    
    
    private StepInputData stepInputData;
    
    private ArrayList<MappedRead> mappedReads ;
    int[][] startPositions = new int[169100][9];
    int[][] readLengths = new int[169100][9];
    
    int[][] mrnaStartPositions = new int[4665][9];
    int[][] mrnaStopPositions  = new int[4665][9];
    
    int[] chrFeatureCount = new int[10];
    
    public StepAnalyzeReadPositions()
    {
       
    }
    
    
    /**
     * read GFF file containing feature information
     * 
     * @param filename
     * @throws IOException 
     */
    public void readGFF(String filename) throws IOException
    {
        String line = null;
        int lineCount = 0;
        BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
	while ((line = br.readLine()) != null) {
            logger.debug(line);

            String chrName = line.split("\t")[0];
            chrName = chrName.replace("chr", "");
            int chr = 0;
            
            if (chrName.equals("M") == false) {
                chr = Integer.parseInt(chrName);
            }
            
            String startPos = line.split("\t")[3];
            String stopPos = line.split("\t")[4];
            
            mrnaStartPositions[chrFeatureCount[chr]][chr] = Integer.parseInt(startPos);
            mrnaStopPositions[chrFeatureCount[chr]][chr] = Integer.parseInt(stopPos);
            chrFeatureCount[chr] ++;
            
            
        }        
    }
    
    
    /**
     * parse out SAM file to retrieve start and stop information for each read
     * 
     * @param filename
     * @throws IOException 
     */
    @Override
    public void execute() 
    {
        try{
            stepInputData.verifyInputData();            
        }
        catch(IOException exIO){
            logger.info("exception parsing InputData" + exIO);
        }

        String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        String samLine = null;
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            
            String positionFile  = posAnalysisFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", positionsExtension);
            String samInputFile = pathToData + FileSeparator + inFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", infileExtension);
            
            logger.info("sam input file is " + samInputFile);
            logger.info("results will be written to " + positionFile);
            
            samLine = null;
            int matchCount5 = 0;
            int matchCount3 = 0;
            int preMatchCount5 = 0;
            int preMatchCount3 = 0;
            int totalCounts = 0;
        
            int lineCount = 0;
            try{
                BufferedReader brSAM = new BufferedReader(new FileReader(new File(samInputFile)));
                while ((samLine = brSAM.readLine()) != null) {
                    
                    logger.debug(samLine);
                    SAMEntry e;
                    SAMEntry samEntry = new SAMEntry(samLine);
                    mappedReads.add(new MappedRead(samEntry.getStartPos(), samEntry.getEndPos(), samEntry.getqName()));

                    if (samLine.startsWith("@")) continue;

                    String qname = samLine.split("\t")[0];

                    int flag = Integer.parseInt(samLine.split("\t")[1]);
                    String rname = samLine.split("\t")[2];
                    if (rname.equals("*")) continue;
                    rname = rname.replace("chr", "");
                    int chr = 0;

                    if (rname.equals("M") == false) {
                        chr = Integer.parseInt(rname);
                    }


                    int pos = Integer.parseInt(samLine.split("\t")[3]);

                    int mapq = Integer.parseInt(samLine.split("\t")[4]);
                    String cigar = samLine.split("\t")[5];
                    String rnext = samLine.split("\t")[6];
                    int pnext = Integer.parseInt(samLine.split("\t")[7]);
                    int tlen = Integer.parseInt(samLine.split("\t")[8]);
                    String sequence = samLine.split("\t")[9];
                    String qual = samLine.split("\t")[10];

                    if ((flag & 0x4) == 0x4) continue;

                    startPositions[lineCount][chr] = pos;
                    readLengths[lineCount][chr] = Integer.parseInt(cigar.replace("M",""));
                    logger.debug(qname + ":\t" + "(" + flag + ")\t" + pos + "\t" + pos + tlen + "\t" + sequence.length());

                    lineCount ++;
                }

            brSAM.close();

            }
            catch(IOException smIO){

            }
        
            try{
                for(int currentChr=0; currentChr<9; currentChr++){
                    String chrInMrnaFileName = samInputFile.replace(".sam", ".chr" + currentChr + ".inreads.tab");
                    String chrOutMrnaFileName = samInputFile.replace(".sam", ".chr" + currentChr + ".outreads.tab");
                    logger.info(chrInMrnaFileName);
                    BufferedWriter bwInReads = new BufferedWriter(new FileWriter(new File(chrInMrnaFileName)));
                    BufferedWriter bwOutReads = new BufferedWriter(new FileWriter(new File(chrOutMrnaFileName)));
                        lineCount = 0;
                        for(int l=0; l<169100; l++){
        //                    logger.debug(startPositions[lineCount][i]);
                            if (startPositions[lineCount][currentChr] != 0) {
                                // does this read hit an mRNA feature?
                                for(int currentMrna=0;currentMrna<chrFeatureCount[currentChr];currentMrna++){
                                    if(startPositions[lineCount][currentChr] > mrnaStartPositions[currentMrna][currentChr] 
                                      && startPositions[lineCount][currentChr] < mrnaStopPositions[currentMrna][currentChr]){
                                        int readStartPosInMrna = startPositions[lineCount][currentChr] - mrnaStartPositions[currentMrna][currentChr];
                                        bwInReads.write(readStartPosInMrna + "\n");
                                    }
                                    else{
                                        String lineOut = startPositions[lineCount][currentChr] + "\t" + readLengths[lineCount][currentChr] + "\n";
                                        bwOutReads.write(lineOut);
                                    }
                                }
        //                        String lineOut = startPositions[lineCount][currentChr] + "\t" + readLengths[lineCount][currentChr] + "\n";
                            }
                            lineCount ++;            
                        }
                    bwOutReads.close();
                    bwInReads.close();
                }
            }
            catch(IOException smIO){

            }
           
        
        }
    }
    
    public static void main(String args[]){
        StepAnalyzeReadPositions readPos = new StepAnalyzeReadPositions();
        try{
            readPos.readGFF("/home/sr/research/sweden/msy2.mrna.gff");
            readPos.execute();
//            readPos.execute("/home/sr/research/sweden/2726_S10.trim.gen.PB1.sam");        
        }
        catch(IOException exIO){
            logger.info("error reading SAM file");
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
