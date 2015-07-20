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
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 *
 * @author sr
 */
public class StepAnalyzeReadPositions {
    
    
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
    
    static Logger logger = LogManager.getLogger();
    
    private StepInputData stepInputData;
    
    //private int[] startPositions;
    //private int[] stopPositions;
    //private int[] readLengths;
    
    int[][] startPositions = new int[169100][9];
    int[][] readLengths = new int[169100][9];
    
    int[][] mrnaStartPositions = new int[4665][9];
    int[][] mrnaStopPositions  = new int[4665][9];
    
    int[] chrFeatureCount = new int[10];
//    int[][] startPositions  = new int[40][9];
//    int[][] readLengths     = new int[40][9];
    
    public StepAnalyzeReadPositions()
    {
       
    }
    
    
    public void readGFF(String filename) throws IOException
    {
        String line = null;
        int lineCount = 0;
        BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
	while ((line = br.readLine()) != null) {
            //logger.debug(line);
            String tokens[] = line.split("\t");

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
    
    
    public void parseSAM(String filename) throws IOException
    {
        String line = null;
        int lineCount = 0;
        BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
	while ((line = br.readLine()) != null) {
            //logger.debug(line);

            if (line.startsWith("@")) continue;

            String qname = line.split("\t")[0];

            int flag = Integer.parseInt(line.split("\t")[1]);
            String rname = line.split("\t")[2];
            if (rname.equals("*")) continue;
            rname = rname.replace("chr", "");
            int chr = 0;
            
            if (rname.equals("M") == false) {
                chr = Integer.parseInt(rname);
            }


            int pos = Integer.parseInt(line.split("\t")[3]);

            int mapq = Integer.parseInt(line.split("\t")[4]);
            String cigar = line.split("\t")[5];
            String rnext = line.split("\t")[6];
            int pnext = Integer.parseInt(line.split("\t")[7]);
            int tlen = Integer.parseInt(line.split("\t")[8]);
            String sequence = line.split("\t")[9];
            String qual = line.split("\t")[10];

            if ((flag & 0x4) == 0x4) continue;

            startPositions[lineCount][chr] = pos;
            readLengths[lineCount][chr] = Integer.parseInt(cigar.replace("M",""));
            //logger.debug(qname + ":\t" + "(" + flag + ")\t" + pos + "\t" + pos + tlen + "\t" + sequence.length());

            lineCount ++;
	}
 
	br.close();
        
        
        
        
        for(int currentChr=0; currentChr<9; currentChr++){
            String chrInMrnaFileName = filename.replace(".sam", ".chr" + currentChr + ".inreads.tab");
            String chrOutMrnaFileName = filename.replace(".sam", ".chr" + currentChr + ".outreads.tab");
            logger.info(chrInMrnaFileName);
            BufferedWriter bwInReads = new BufferedWriter(new FileWriter(new File(chrInMrnaFileName)));
            BufferedWriter bwOutReads = new BufferedWriter(new FileWriter(new File(chrOutMrnaFileName)));
                lineCount = 0;
                for(int l=0; l<169100; l++){
//                    logger.info(startPositions[lineCount][i]);
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
    
    public static void main(String args[]){
        StepAnalyzeReadPositions readPos = new StepAnalyzeReadPositions();
        try{
            readPos.readGFF("/home/sr/research/sweden/msy2.mrna.gff");
            readPos.parseSAM("/home/sr/research/sweden/2726_S10.trim.gen.PB1.sam");        
        }
        catch(IOException exIO){
            logger.info("error reading SAM file");
        }
           
    }
      
}
