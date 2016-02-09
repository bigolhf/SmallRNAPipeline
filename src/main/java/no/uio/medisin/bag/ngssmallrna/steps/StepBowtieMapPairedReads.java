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
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;


/**
 *  Map paired reads to contaminants, RNA, Reference Genome and output in SAM format
 * 
 *   Input is a collapsed FASTA file
 * 
 * @author sr
 */

public class StepBowtieMapPairedReads extends NGSStep{
    
    static Logger                       logger                      = LogManager.getLogger();
    static  String                      FileSeparator               = System.getProperty("file.separator");
    
    private static final String         inFolder                    = "collapsed";
    private static final String         abundReadsOutFolder         = "bowtie_abundant_mapped";
    private static final String         genomeReadsOutFolder        = "bowtie_genome_mapped";
    
    
    private static final String         infileExtension             = ".trim.clp.fasta";
    private static final String         fastqAbundantAlnExtension   = ".trim.clp.abun.fasta";
    private static final String         fastqAbundantUnAlnExtension = ".trim.clp.notabun.fasta";
    private static final String         samAbundantAlnExtension     = ".trim.clp.abun.sam";
    private static final String         fastqGenomeAlnExtension     = ".trim.clp.gen.fasta";
    private static final String         fastqGenomeUnAlnExtension   = ".trim.clp.unmap.fasta";
    private static final String         samGenomeAlnExtension       = ".trim.clp.gen.sam";
    private static final String         mappingResultExtension      = ".trim.clp.gen.mapping.txt";
    
    
    
    private StepInputData               stepInputData;
    private StepResultData              stepResultData;
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepBowtieMapPairedReads(StepInputData sid){
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
        
            bowtieMapReadsParams.put("bowtieMapGenomeRootFolder", this.getGenomeRootFolder());
            bowtieMapReadsParams.put("bowtieMapAlignMode", this.getBowtieMappingAlignMode());
            bowtieMapReadsParams.put("bowtieMapNoOfMismatches", this.getBowtieMappingNoOfMismatches());
            bowtieMapReadsParams.put("bowtieNoOfThreads", this.getBowtieMappingNoOfThreads());
            bowtieMapReadsParams.put("bowtieReferenceGenome", this.getBowtieMappingReferenceGenome());

            bowtie <ref_genome> <fasta_file>
        */
        String mappingCmd = (String) stepInputData.getStepParams().get("bowtieMappingCommand");
        logger.info("Mapping command is " + mappingCmd);
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            ArrayList<String> cmd = new ArrayList<>();
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                String pathToData = stepInputData.getProjectRoot() + FileSeparator + stepInputData.getProjectID();
                
                /*
                bowtie -a -v 2 e_coli --suppress 1,5,6,7 -c ATGCATCATGCGCCA

                use:
                    -f input files are FASTA format

                    -v option (which ignores quality values) since we are using FASTA files,
                    --best (to order the matches)
                    -m 2 because we only want reads that map to a unique location 
                        (we allow 2, because some miRNAs have two locations)

                    --al aligned reads
                    --un unaligned reads
                    --sam SAM file name

                */      
                
                /*
                    Map Abundant Reads
                */
                    String fastqRaw2 = pathToData + FileSeparator + "fastq_files" + FileSeparator + sampleData.getFastqFile1();                    
                cmd.add(mappingCmd);
                String pathToBowtieIndex = stepInputData.getStepParams().get("bowtieMapGenomeRootFolder") 
                    + FileSeparator + stepInputData.getStepParams().get("bowtieReferenceGenome") + "/Sequence/AbundantSequences/abundant";
                cmd.add(pathToBowtieIndex);
                
                String fastqTrimmedInputFile = pathToData + FileSeparator + inFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", infileExtension);
                fastqTrimmedInputFile = fastqTrimmedInputFile.replace(FileSeparator + FileSeparator, FileSeparator).trim();
                cmd.add("-f");
                cmd.add(fastqTrimmedInputFile);
                
                cmd.add("-v" + stepInputData.getStepParams().get("bowtieMapNoOfMismatches"));
                cmd.add("--best");
                cmd.add("-m 2");
                

                String AbundantMapOutputFolder = pathToData + FileSeparator + abundReadsOutFolder;
                AbundantMapOutputFolder = AbundantMapOutputFolder.replace(FileSeparator + FileSeparator, FileSeparator).trim();
                Boolean fA = new File(AbundantMapOutputFolder).mkdir();       
                if (fA) logger.info("created output folder <" + AbundantMapOutputFolder + "> for results" );

                String fastqAbundantAln     = pathToData + FileSeparator + abundReadsOutFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", fastqAbundantAlnExtension);
                String fastqAbundantUnAln   = pathToData + FileSeparator + abundReadsOutFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", fastqAbundantUnAlnExtension);
                String samAbundantAln       = pathToData + FileSeparator + abundReadsOutFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", samAbundantAlnExtension);
                cmd.add("--al " + fastqAbundantAln);
                cmd.add("--un " + fastqAbundantUnAln);
                cmd.add("--sam "  + samAbundantAln);
                cmd.add("");
                

                String cmdBowtieMapAbunReads = StringUtils.join(cmd, " ");
                cmdBowtieMapAbunReads = cmdBowtieMapAbunReads.replace(FileSeparator + FileSeparator, FileSeparator);
                logger.info("Bowtie Map Abundant Reads command:\t" + cmdBowtieMapAbunReads);

                Runtime rt = Runtime.getRuntime();
                Process proc = rt.exec(cmdBowtieMapAbunReads);
                BufferedReader brAStdin  = new BufferedReader(new InputStreamReader(proc.getInputStream()));
                BufferedReader brAStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
                
                    String line = null;
                    logger.info("<OUTPUT>");
                    while ((line = brAStdin.readLine()) != null)
                        logger.info(line);
                    logger.info("</OUTPUT>");

                    logger.info("<ERROR>");
                    int skipCount = 0;
                    ArrayList<String > mapAbunStdErr = new ArrayList<>();
                    while ( (line = brAStdErr.readLine()) != null){
                        if(line.contains("Warning: Skipping") && line.contains("less than")){
                            skipCount++;
                        }else{
                            logger.info(line);
                            mapAbunStdErr.add(line);
                        }
                    }
                    // need to parse the output from Bowtie to get the mapping summary
                    logger.info(skipCount + " lines were skipped because the read was too short");
                    logger.info("</ERROR>");
                
                    int exitVal = proc.waitFor();            
                    logger.info("Process exitValue: " + exitVal);   
                
                brAStdin.close();
                brAStdErr.close();
                
                
                
                /*
                    Map Genome Reads
                */
                cmd = new ArrayList<>();
                cmd.add(mappingCmd);
                String pathToBowtieGenomeIndex = stepInputData.getStepParams().get("bowtieMapGenomeRootFolder") 
                    + FileSeparator + stepInputData.getStepParams().get("bowtieReferenceGenome") + "/Sequence/BowtieIndex/genome";
                cmd.add(pathToBowtieGenomeIndex);
                
//                fastqTrimmedInputFile = fastqAbundantUnAln;
                cmd.add("-f");
                cmd.add(fastqAbundantUnAln);
                
                cmd.add("-v" + stepInputData.getStepParams().get("bowtieMapNoOfMismatches"));
                cmd.add("--best");
                cmd.add("-m 2");
                

                String genomeMapOutputFolder = pathToData + FileSeparator + genomeReadsOutFolder;
                genomeMapOutputFolder = genomeMapOutputFolder.replace(FileSeparator + FileSeparator, FileSeparator).trim();
                Boolean fG = new File(genomeMapOutputFolder).mkdir();       
                if (fG) logger.info("created genome map output folder <" + genomeMapOutputFolder + "> for results" );

                String fastqGenomeAln     = genomeMapOutputFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", fastqGenomeAlnExtension);
                String fastqGenomeUnAln   = genomeMapOutputFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", fastqGenomeUnAlnExtension);
                String samGenomeAln       = genomeMapOutputFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", samGenomeAlnExtension);
                cmd.add("--al " + fastqGenomeAln);
                cmd.add("--un " + fastqGenomeUnAln);
                cmd.add("--sam "  + samGenomeAln);
                

                String cmdBowtieMapGenomeReads = StringUtils.join(cmd, " ");
                cmdBowtieMapGenomeReads = cmdBowtieMapGenomeReads.replace(FileSeparator + FileSeparator, FileSeparator);
                logger.info("Bowtie Map Genome Reads command:\t" + cmdBowtieMapGenomeReads);

                Runtime rtGenMap = Runtime.getRuntime();
                Process procGenMap = rtGenMap.exec(cmdBowtieMapGenomeReads);
                BufferedReader brGStdin  = new BufferedReader(new InputStreamReader(procGenMap.getInputStream()));
                BufferedReader brGStdErr = new BufferedReader(new InputStreamReader(procGenMap.getErrorStream()));
                
                    String gLine = null;
                    logger.info("<OUTPUT>");
                    while ( (gLine = brGStdin.readLine()) != null)
                        logger.info(line);
                    logger.info("</OUTPUT>");

                    logger.info("<ERROR>");
                    skipCount = 0;
                    ArrayList<String > mapGenStdErr = new ArrayList<>();
                    while ( (line = brGStdErr.readLine()) != null){
                        if(line.contains("Warning: Skipping") && line.contains("less than")){
                            skipCount++;
                        }else{
                            logger.info(line);
                            mapGenStdErr.add(line);
                        }
                    }
                    // need to parse the output from Bowtie to get the mapping summary
                    logger.info(skipCount + " lines were skipped because the read was too short");
                    logger.info("</ERROR>");
                
                    int gExitVal = procGenMap.waitFor();            
                    logger.info("Process exitValue: " + gExitVal);   
                
                brAStdin.close();
                brAStdErr.close();
                

                String mappingOutputFile    = genomeMapOutputFolder + FileSeparator + sampleData.getFastqFile1().replace(".fastq", mappingResultExtension);
                BufferedWriter bwMO = new BufferedWriter(new FileWriter(new File(mappingOutputFile)));
                    for(String mapLine: mapGenStdErr){
                        bwMO.write(mapLine + "\n");
                    }
                    bwMO.write("\n\n" + "+" + StringUtils.repeat("-", 60) + "+" + "\n");
                    
                    String fastqRaw = pathToData + FileSeparator + "fastq_files" + FileSeparator + sampleData.getFastqFile1();                    
                    int rawReadsIn = 0;
                    BufferedReader brFQ = new BufferedReader(new FileReader(fastqRaw));
                        while (brFQ.readLine() != null) rawReadsIn++;
                    brFQ.close();

                bwMO.write(fastqRaw + "\n");
                bwMO.write(fastqTrimmedInputFile + "\n");
                bwMO.write(fastqGenomeAln + "\n");
                bwMO.write(fastqAbundantAln + "\n");
                bwMO.write(fastqGenomeUnAln + "\n");

                    // Input 
                    int totalInputReads = 0;
                    String fqLine = "";
                    BufferedReader brIR = new BufferedReader(new FileReader(new File(fastqTrimmedInputFile)));
                        while((fqLine = brIR.readLine())!=null){
                            totalInputReads += Integer.parseInt(fqLine.substring(1).split("-")[1]);
                            
                            brIR.readLine();
                        }
                    brIR.close();
                    bwMO.write("total input reads = " + rawReadsIn + "\n");
                    
                    // Mapped
                    int totalMappedReads = 0;
                    fqLine = "";
                    BufferedReader brMR = new BufferedReader(new FileReader(new File(fastqGenomeAln)));
                        while((fqLine = brMR.readLine())!=null){
                            totalMappedReads += Integer.parseInt(fqLine.substring(1).split("-")[1]);
                            brMR.readLine();
                        }
                    brMR.close();
                    bwMO.write("total mapped reads = " + totalMappedReads + "\n");
                    
                    // Abundant
                    int totalAbundantReads = 0;
                    BufferedReader brAR = new BufferedReader(new FileReader(new File(fastqAbundantAln)));
                        while((fqLine = brAR.readLine())!=null){
                            totalAbundantReads += Integer.parseInt(fqLine.substring(1).split("-")[1]);
                            brAR.readLine();
                        }
                    brAR.close();
                    bwMO.write("total abundant reads = " + totalAbundantReads  + "\n");

                    // Unmapped
                    int totalUnmappedReads = 0;
                    BufferedReader brUR = new BufferedReader(new FileReader(new File(fastqGenomeUnAln)));
                        while((fqLine = brUR.readLine())!=null){
                            totalUnmappedReads += Integer.parseInt(fqLine.substring(1).split("-")[1]);
                            brUR.readLine();
                        }
                    brUR.close();
                    bwMO.write("total unmapped reads = " + totalUnmappedReads  + "\n");
                    
                    bwMO.write("length filtered reads = " 
                            + (rawReadsIn - totalMappedReads - totalAbundantReads - totalUnmappedReads));

                    bwMO.write("\n\n" + "+" + StringUtils.repeat("-", 60) + "+" + "\n");
                    
                bwMO.close();
            }
            catch(IOException|InterruptedException ex){
                logger.error("error executing Bowtie Mapping command\n");
                logger.error(cmd);
                logger.error(ex.toString());
            }
        }
        
        
    }
    
    
            
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
    public void verifyOutputData(){
        
    }
}
