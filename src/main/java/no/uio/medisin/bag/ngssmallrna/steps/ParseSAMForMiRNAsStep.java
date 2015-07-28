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
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import no.uio.medisin.bag.ngssmallrna.pipeline.MiRNAFeature;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;


/**
 *  parse SAM file to extract and process the miRNA reads
 * 
 *   Input is a SAM file
 * 
 * @author sr
 */

public class ParseSAMForMiRNAsStep extends NGSStep{
    
    static Logger                       logger                      = LogManager.getLogger();
    static  String                      FileSeparator               = System.getProperty("file.separator");
    
    private static final String         inFolder                    = "bowtie_genome_mapped";
    private static final String         abundReadsOutFolder         = "bowtie_abundant_mapped";
    private static final String         genomeReadsOutFolder        = "bowtie_genome_mapped";
    
    
    private static final String         infileExtension             = ".trim.clp.gen.sam";
    private static final String         fastqAbundantAlnExtension   = ".trim.clp.abun.fasta";
    private static final String         fastqAbundantUnAlnExtension = ".trim.clp.notabun.fasta";
    private static final String         samAbundantAlnExtension     = ".trim.clp.abun.sam";
    private static final String         fastqGenomeAlnExtension     = ".trim.clp.gen.sam";
    private static final String         fastqGenomeUnAlnExtension   = ".trim.clp.unmap.fasta";
    private static final String         samGenomeAlnExtension       = ".trim.clp.gen.sam";
    
    
    
    private StepInputData               stepInputData;
    private StepResultData              stepResultData;
    

    private List<MiRNAFeature>          miRNAList                   = new ArrayList<>();
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public ParseSAMForMiRNAsStep(StepInputData sid){
        try{
            stepInputData = sid;
            stepInputData.verifyInputData();
            
        }
        catch(IOException exIO){
            
        }
    }
    
    @Override
    public void execute(){
        /*
            parseSAMmiRNAsParams.put("bleed", this.getSamParseForMiRNAsBleed());
            parseSAMmiRNAsParams.put("miRBaseHostGFFFile", this.getMiRBaseHostGFF());
        
        */
        
        try{
            this.loadMiRBaseData((String) stepInputData.getStepParams().get("miRBaseHostGFFFile"));
        }
        catch(IOException ex){
            logger.error("error reading miRBase reference file <" + (String) stepInputData.getStepParams().get("miRBaseHostGFFFile") + ">\n" + ex.toString());
        }
    
        
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                int bleed = (int) stepInputData.getStepParams().get("bleed");
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                
                String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();                
                String samInputFile = pathToData + FileSeparator + inFolder + FileSeparator + sampleData.getDataFile().replace(".fastq", infileExtension);
                int matchCount5 = 0;
                int matchCount3 = 0;
                int preMatchCount5 = 0;
                int preMatchCount3 = 0;
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
                            
                            if (this.doesReadOverlapKnownMiRNA(startPos, endPos, chr, strand, bleed) != null){
                                List<String> mutations = new ArrayList<>();
                                Matcher match = Pattern.compile("[0-9]+|[a-z]+|[A-Z]+").matcher(samLine.split("\t")[12].split(":")[2]);
                                String outputString = samLine.split("\t")[0] + ":" + samLine.split("\t")[12] + ": ";
                                while (match.find()) {
                                    mutations.add(match.group());
                                    outputString = outputString.concat(match.group() + "|");
                                }
                                logger.info(outputString);
                                
                            }
                            
                            
                            
                            
                            
                        }
                    }
                    logger.info((matchCount5 + matchCount3) + " reads (" + matchCount5 + " 5'" + "/" + matchCount3 + " 3' ) were mapped");
                brSAM.close();
                
                
            }
            catch(IOException ex){
                logger.error("error executing AdapterTrimming command\n" + ex.toString());
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
        
        for(MiRNAFeature miRBaseEntry: this.miRNAList){
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
     * load miRNA specs (name and Chromosome position) from GFF file
     * downloaded from miRBase
     * Because different releases of miRBase use different releases of
     * reference genome, we have to track both miRBase and genome reference IDs
     * 
     * @param miRBaseGFFFile    String : absolute path to file
     * @throws IOException
     * 
     */
    public void loadMiRBaseData(String miRBaseGFFFile) throws IOException{
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
                this.miRNAList.add(new MiRNAFeature(name, chr, startPos, endPos, strand, id, parent));
            }
        brMiR.close();
        logger.info("read " + miRNAList + "miRNA entries");
    }
    
    
    
    @Override
    public void outputResultData(){
        
    }
}
