/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.jmirpara.SimpleSeq;
import no.uio.medisin.bag.ngssmallrna.pipeline.ReferenceDataLocations;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;




/**
 *  Adapter Trimming Step
 *  perform QC of FASTQ file list.
 * 
 *   Input is a FASTQ file
 *   Output is a quality report
 * 
 * 
 * @author sr
 */

public class StepBSMap extends NGSStep implements NGSBase{
    
    private static Logger               logger                      = LogManager.getLogger();
    
    private static final String         STEP_ID_STRING              = "BSMapReads";

    private static final String         ID_SOFTWARE                 = "pathToBSmap:";
    private static final String         ID_REF_GENOME               = "host:";
    private static final String         ID_MISMATCHES               = "noOfMismatches:";
    private static final String         ID_SEED_SIZE                = "seedSize:";
    //private static final String         ID_MISMATCH                 = "mismatches";
    private static final String         ID_THREEMAPPING             = "threeNucleoideMapping";
    private static final String         ID_GAPSIZE                  = "gapSize";
    private static final String         ID_KMERCUTOFF               = "kmerCutoffRatio:";
    private static final String         ID_REPORT_REPEATS           = "reportRepeatsValue:";
    private static final String         ID_QUAL_THRESHOLD_VALUE     = "qualityThresholdTrimValue:";
    private static final String         ID_LOW_QUAL_VALUE           = "lowQualityFilterValue:";
    private static final String         ID_THREADS                  = "noOfThreads:";
    private static final String         ID_MAX_INSERT_SIZE          = "maxInsertSize:";
    private static final String         ID_MIN_INSERT_SIZE          = "minInsertSize:";
    private static final String         ID_MAP_FIRST_N_NUCS         = "mapFirstNucleotides:";
    private static final String         ID_GENOME_INDEX_INTERVAL    = "genomeIndexInterval:";
    private static final String         ID_ADAPTERSEQ               = "adapterSequence:";
    private static final String         ID_TRIM_ADAPTERSEQ          = "trimAdapterSequence:";
    private static final String         ID_INCLUDE_REF_SEQ          = "includeRefSeq:";
    private static final String         ID_SKIP_SAM_HEADER          = "skipSAMHeader:";
    private static final String         ID_REPORT_UNMAPPED_READS    = "reportUnmappedReads:";
    private static final String         ID_START_AT_READ_N          = "startAtThisRead:";
    private static final String         ID_END_AT_READ_N            = "endtAtThisRead:";
    private static final String         ID_DIGESTIONSITE            = "digestionSite:";
    private static final String         ID_RANDOMSEED               = "randomSeed:";
    private static final String         ID_MAPPINGSTRAND            = "mappingStrand:";
    private static final String         ID_TRANSITIONMAP            = "transitionMap:";
    private static final String         ID_MESSAGELEVEL             = "messageLevel:";
    
    
    private static final String         INFILE_EXTENSION            = ".fastq";
    private static final String         OUTFILE_EXTENSION           = ".fastq";

    /**
     * @return the logger
     */
    public static Logger getLogger() {
        return logger;
    }

    /**
     * @param aLogger the logger to set
     */
    public static void setLogger(Logger aLogger) {
        logger = aLogger;
    }

    
    private String                      pathToSoftware              = "";

    private             String          pathToBSMap                 = "";
    private             int             seedSize                    = 0;
    private             String          AlignMode                   = "";
    private             double          noOfMismatches              = 0.08;
    private             int             gapSize                     = 0;
    private             Boolean         threeNucMapping             = false;
    private             int             noOfThreads                 = 4;
    private             String          rootDataFolder              = "";
    private             String          ReferenceGenome             = "";
    private             double          kmerCutoffRatio             = 1e-06;
    private             int             reportRepeatsValue          = 1;
    private             int             qualityThresholdTrimValue   = 0;
    private             int             lowQualityFilterValue       = 5;
    private             int             maxInsertSize               = 500;
    private             int             minInsertSize               = 28;
    private             int             mapFirstNucleotides         = 0;
    private             int             genomeIndexInterval         = 4;
    private     ArrayList<String>       adapterSequences;
    private             Boolean         trimAdapterSequence         = false;
    private             Boolean         includeRefSeq               = false;
    private             Boolean         skipSAMHeader               = false;
    private             Boolean         reportUnmappedReads         = false;
    private             int             startAtThisRead             = 1;
    private             int             endtAtThisRead              = -1;
    private             String          digestionSite               = "";
    private             int             randomSeed                  = 0;
    private             int             mappingStrand               = 0;
    private             String          transitionMap               = "TC";
    private             int             messageLevel                = 1;
     

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepBSMap(StepInputData sid){
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

        getLogger().info(STEP_ID_STRING + ": verify configuration data");
        
        /*
            the following parameters have to be specfied in the configuration file for the 
            analysis to proceed
        */
        if(configData.get(ID_SOFTWARE)==null) {
            getLogger().info("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
            getLogger().error("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            getLogger().info("<" + configData.get(ID_REF_GENOME) + "> : Missing Definition in Configuration File");
            getLogger().error("<" + configData.get(ID_REF_GENOME) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_REF_GENOME) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_THREADS)==null) {
            getLogger().info("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
            getLogger().error("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_DIGESTIONSITE)==null) {
            getLogger().info("<" + configData.get(ID_DIGESTIONSITE) + "> : Missing Definition in Configuration File");
            getLogger().error("<" + configData.get(ID_DIGESTIONSITE) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_DIGESTIONSITE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_TRIM_ADAPTERSEQ)==null) {
            getLogger().info("<" + configData.get(ID_TRIM_ADAPTERSEQ) + "> : Missing Definition in Configuration File");
            getLogger().error("<" + configData.get(ID_TRIM_ADAPTERSEQ) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_TRIM_ADAPTERSEQ) + "> : Missing Definition in Configuration File");
        }
        


        String chk;

        chk = checkParameter("Integer", ID_SEED_SIZE, (String)configData.get(ID_SEED_SIZE), "8", "16", logger);
        if(chk!=null)
            this.setSeedSize((Integer)configData.get(ID_SEED_SIZE));
            
        chk = checkParameter("Double", ID_MISMATCHES, (String)configData.get(ID_MISMATCHES), "0", "15", logger);
        if(chk!=null)
            this.setNoOfMismatches((Double)configData.get(ID_MISMATCHES));
            
        chk = checkParameter("Integer", ID_GAPSIZE, (String)configData.get(ID_GAPSIZE), "0", "NA", logger);
        if(chk!=null)
            this.setGapSize((Integer)configData.get(ID_GAPSIZE));
                        
        chk = checkParameter("Boolean", ID_THREEMAPPING, (String)configData.get(ID_THREEMAPPING), "0", "NA", logger);
        if(chk!=null)
            this.setThreeNucMapping((Boolean)configData.get(ID_THREEMAPPING));
                        
        chk = checkParameter("Integer", ID_THREADS, (String)configData.get(ID_THREADS), "0", "12", logger);
        if(chk!=null)
            this.setNoOfThreads((Integer)configData.get(ID_THREADS));
            
        chk = checkParameter("Double", ID_KMERCUTOFF, (String)configData.get(ID_KMERCUTOFF), "0", "12", logger);
        if(chk!=null)
            this.setKmerCutoffRatio((Double)configData.get(ID_KMERCUTOFF));
            
        chk = checkParameter("Integer", ID_REPORT_REPEATS, (String)configData.get(ID_REPORT_REPEATS), "0", "2", logger);
        if(chk!=null)
            this.setReportRepeatsValue((Integer)configData.get(ID_REPORT_REPEATS));
            
        chk = checkParameter("Integer", ID_QUAL_THRESHOLD_VALUE, (String)configData.get(ID_QUAL_THRESHOLD_VALUE), "0", "40", logger);
        if(chk!=null)
            this.setQualityThresholdTrimValue((Integer)configData.get(ID_QUAL_THRESHOLD_VALUE));
            
        chk = checkParameter("Integer", ID_LOW_QUAL_VALUE, (String)configData.get(ID_LOW_QUAL_VALUE), "0", "NA", logger);
        if(chk!=null)
            this.setLowQualityFilterValue((Integer)configData.get(ID_LOW_QUAL_VALUE));
            
        chk = checkParameter("Integer", ID_MAX_INSERT_SIZE, (String)configData.get(ID_MAX_INSERT_SIZE), "0", "NA", logger);
        if(chk!=null)
            this.setMaxInsertSize((Integer)configData.get(ID_MAX_INSERT_SIZE));
            
        chk = checkParameter("Integer", ID_MIN_INSERT_SIZE, (String)configData.get(ID_MIN_INSERT_SIZE), "0", "NA", logger);
        if(chk!=null)
            this.setMinInsertSize((Integer)configData.get(ID_MIN_INSERT_SIZE));
            
        chk = checkParameter("Integer", ID_MAP_FIRST_N_NUCS, (String)configData.get(ID_MAP_FIRST_N_NUCS), "0", "NA", logger);
        if(chk!=null)
            this.setMapFirstNucleotides((Integer)configData.get(ID_MAP_FIRST_N_NUCS));
            
        chk = checkParameter("Integer", ID_GENOME_INDEX_INTERVAL, (String)configData.get(ID_GENOME_INDEX_INTERVAL), "0", "NA", logger);
        if(chk!=null)
            this.setGenomeIndexInterval((Integer)configData.get(ID_GENOME_INDEX_INTERVAL));
 
        chk = checkParameter("Boolean", ID_TRIM_ADAPTERSEQ, (String)configData.get(ID_TRIM_ADAPTERSEQ), "NA", "NA", logger);
        if(chk!=null && this.getTrimAdapterSequence()==true && this.checkAdaptorSequences()!=null){
            this.setAdapterSequences((ArrayList<String>)configData.get(ID_TRIM_ADAPTERSEQ));
        }

        chk = checkParameter("Boolean", ID_INCLUDE_REF_SEQ, (String)configData.get(ID_INCLUDE_REF_SEQ), "NA", "NA", logger);
        if(chk!=null)
            this.setIncludeRefSeq((Boolean)configData.get(ID_INCLUDE_REF_SEQ));
        
        chk = checkParameter("Boolean", ID_SKIP_SAM_HEADER, (String)configData.get(ID_SKIP_SAM_HEADER), "NA", "NA", logger);
        if(chk!=null)
            this.setSkipSAMHeader((Boolean)configData.get(ID_SKIP_SAM_HEADER));
        
        chk = checkParameter("Boolean", ID_REPORT_UNMAPPED_READS, (String)configData.get(ID_REPORT_UNMAPPED_READS), "NA", "NA", logger);
        if(chk!=null)
            this.setReportUnmappedReads((Boolean)configData.get(ID_REPORT_UNMAPPED_READS));       
        
        chk = checkParameter("Integer", ID_START_AT_READ_N, (String)configData.get(ID_START_AT_READ_N), "0", "NA", logger);
        if(chk!=null)
            this.setStartAtThisRead((Integer)configData.get(ID_START_AT_READ_N));
            
        chk = checkParameter("Integer", ID_END_AT_READ_N, (String)configData.get(ID_END_AT_READ_N), "0", "NA", logger);
        if(chk!=null)
            this.setEndtAtThisRead((Integer)configData.get(ID_END_AT_READ_N));
        
        if(this.checkDigestionSite((String)configData.get(ID_DIGESTIONSITE))){
            this.setDigestionSite((String)configData.get(ID_DIGESTIONSITE));
        }

        checkParameter("Integer", ID_RANDOMSEED, (String)configData.get(ID_RANDOMSEED), "0", "NA", logger);
        if(chk!=null)
            this.setRandomSeed((Integer)configData.get(ID_RANDOMSEED));
            
        checkParameter("Integer", ID_MAPPINGSTRAND, (String)configData.get(ID_MAPPINGSTRAND), "0", "NA", logger);
        if(chk!=null)
            this.setMappingStrand((Integer)configData.get(ID_MAPPINGSTRAND));
            
        if(this.checkTransitionMap((String)configData.get(ID_TRANSITIONMAP))){
            this.setTransitionMap((String)configData.get(ID_TRANSITIONMAP));
        }

        checkParameter("Integer", ID_MESSAGELEVEL, (String)configData.get(ID_MESSAGELEVEL), "0", "NA", logger);
        if(chk!=null)
            this.setMessageLevel((Integer)configData.get(ID_MESSAGELEVEL));
            

        getLogger().info("passed");
    }
    
    
    
    /**
     * check the adaptor sequence is comprised of standard DNA nucleotides
     * 
     * @param seq
     * @return 
     */
    private Boolean checkAdaptorSequences() throws Exception{
        
        for(String adapterSeq: this.getAdapterSequences()){
            for (char ch: adapterSeq.toCharArray()) {
                if("acgtACGT".indexOf(ch) == -1){
                    logger.info("Adaptor Sequence <" + adapterSeq + ">" + " can only contain standard nucleotides (ACGT)");            
                    logger.error("Adaptor Sequence <" + adapterSeq + ">" + " can only contain standard nucleotides (ACGT)");   
                    throw new Exception("Adaptor Sequence <" + adapterSeq + ">" + " can only contain standard nucleotides (ACGT)");

                }
            }
            
        }
        return true;
    }
    



    /**
     * Digestion site should contain standard DNA nucleotides and a '-'
     * to mark the digestion location
     * 
     * @param seq
     * @return 
     */
    private Boolean checkDigestionSite(String seq) throws Exception{
        if(seq.contains("-")==false) return false;
        
        String rawSeq = "";
        for (char ch: seq.toCharArray()) {
            if("acgtACGT-".indexOf(ch) == -1){
                logger.info("Digestion Site <" + seq + ">" + " can only contain standard nucleotides (ACGT)");            
                logger.error("Digestion Site <" + seq + ">" + " can only contain standard nucleotides (ACGT)");            
                throw new Exception("Digestion Site <" + seq + ">" + " can only contain standard nucleotides (ACGT)");
            }
            rawSeq = rawSeq + ch;
        }
        
        if(SimpleSeq.complement(seq).equals(seq)==false){
            logger.info("Digestion Site <" + seq + ">" 
                    + " must be palindromic: Complement<" + seq + ">=" + SimpleSeq.complement(seq));            
            logger.error("Digestion Site <" + seq + ">" 
                    + " must be palindromic: Complement<" + seq + ">=" + SimpleSeq.complement(seq));            
            throw new Exception("Digestion Site <" + seq + ">" 
                    + " must be palindromic: Complement<" + seq + ">=" + SimpleSeq.complement(seq));            
        }
        
        return true;
    }
    



    /**
     * Transition Map should comprise a DNA nucleotide pair
     * 
     * @param seq
     * @return 
     */
    private Boolean checkTransitionMap(String seq) throws Exception{
        if(seq.length()!=2) {
            logger.info("Transition Map can only be two nucleotides, <" + seq + "> is "  + seq.length() + " nt");            
            logger.error("Transition Map can only be two nucleotides, <" + seq + "> is "  + seq.length() + " nt");      
            throw new Exception("Transition Map can only be two nucleotides, <" + seq + "> is "  + seq.length() + " nt");
        }
        for (char ch: seq.toCharArray()) {
            if("acgtACGT".indexOf(ch) == -1) {
                logger.info("Transition Map <" + seq + ">" + " can only contain standard nucleotides (ACGT)");            
                logger.error("Transition Map <" + seq + ">" + " can only contain standard nucleotides (ACGT)");            
                throw new Exception("Transition Map <" + seq + ">" + " can only contain standard nucleotides (ACGT)");
            }
        }
        
        return true;
    }
    



    
    /**
     * BSMap single or paired end reads
     * Still need to add ability to include multiple adapter sequences
     * 
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException{
        /* example call.:
            ./bsmap \
            -a ../../../MetylationDataEPIRA/PilotDNAextractionRnBeads/160125_M01334/1-RA5111B-4nQ_S1_L001_R1_001.paired.fastq.gz \
            -b ../../../MetylationDataEPIRA/PilotDNAextractionRnBeads/160125_M01334/1-RA5111B-4nQ_S1_L001_R2_001.paired.fastq.gz \
            -d ../../../../../../Volumes/LaCie_kari/ref/Homo_sapiens_assembly19.fasta \
            -D C-CGG \
            -p 2 \
            -w 100 \
            -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG \
            -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
            -v 2 \
            -n 1 \
            -o testOut29febSAMformat.sam &
        */
        getLogger().info(STEP_ID_STRING + ": execute step");        
        
        String fastqFile1in = "";
        String fastqFile1out = "";
        String fastqFile2in = "";
        String fastqFile2out = "";
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();                
                
                String cmdBSMap = "";   
                ArrayList<String> cmd = new ArrayList<>();
                cmd.add(this.getPathToSoftware());
                
                fastqFile1in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION));
                cmd.add("-a " + fastqFile1in);
                if (sampleData.getFastqFile2() != null){
                    fastqFile2in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(".fastq", INFILE_EXTENSION));
                    cmd.add("-b " + fastqFile2in);
                }
                
                String pathToGenomeFA = this.cleanPath(stepInputData.getDataLocations().getGenomeRootFolder()
                        + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_WHOLE_GENSEQ_FA);
                cmd.add("-d " + pathToGenomeFA);
                    
                /**
                 * not sure about the best way to handle this.
                 * It seems that the files that come of the NSC sequences have a general format
                 * SAMPLENAME_SEQUENCERDATA1_SEQUENCEDATA_..., so we split at the first "_"
                 * 
                 */
                String samAlnFile = this.cleanPath(outFolder + FILESEPARATOR 
                        + sampleData.getFastqFile1().replace(".fastq", "sam").split("_")[0].trim());
                cmd.add("-o " + samAlnFile);
                
                cmd.add("-s " + this.getSeedSize());
                
                /*  
                  if this value is between 0 and 1, it's interpreted as the mismatch rate w.r.t to the read length.
                  otherwise it's interpreted as the maximum number of mismatches allowed on a read, <=15.
                  example: 
                    -v 5 (max #mismatches = 5), 
                    -v 0.1 (max #mismatches = read_length * 10%)
                    default=0.08.
                */
                cmd.add("-v " + this.getNoOfMismatches());
                
                /*
                  using 3-nucleotide mapping approach. (default: off)
                */
                if(this.getThreeNucMapping())
                    cmd.add("-3");
                
                        
                /*
                  gap size, BSMAP only allow 1 continuous gap (insertion or deletion) 
                  with up to 3 nucleotides. default=0
                  gaps will not be allowed within 6nt of the read edges.
                  the number of mismatches of gapped algnment is calculated as #gap_size+#mismatches+1
                */
                cmd.add("-g " + this.getGapSize());
                
                /*
                  set the cut-off ratio for over-represented kmers, default=1e-06
                  e.g.: -k 1e-6 means the top 0.0001% over-represented kmer will be skipped in alignment
                */
                cmd.add("-k " + this.getKmerCutoffRatio());
                
                /*
                  [0,1,2] how to report repeat hits, 0=none(unique hit/pair only); 
                  1=random one; 2=all(large output file size), default=1.
                */                
                cmd.add("-r " + this.getReportRepeatsValue());

                /*
                  quality threshold in trimming 3'end of reads, 0-40, default=0. (no trim)
                */
                cmd.add("-q " + this.getQualityThresholdTrimValue());
                
                /*
                  filter low-quality reads containing >n Ns, default=5
                */
                cmd.add("-f " + this.getLowQualityFilterValue());
                
                cmd.add("-p " + this.getNoOfThreads());
                cmd.add("-x " + this.getMaxInsertSize());
                cmd.add("-m " + this.getMinInsertSize());
                
                /*
                  map the first N nucleotide of the read, default: 0 (map the whole read).
                */
                cmd.add("-L " + this.getMapFirstNucleotides());

                /*
                  index interval (1~16), meaning the reference genome will be indexed every Nbp, default=4. (WGBS mode)
                  For RRBS mode, index_interval is fixed to 1bp and this command line option is neglected.
                 larger index interval uses memory, and slightly reduces mapping sensitivity. (~0.5% difference)
                 for human genome, -I 16 uses ~5GB, compared with ~9GB at the default -I 4.                
                */
                cmd.add("-I " + this.getGenomeIndexInterval());
                
                if(this.getTrimAdapterSequence()){
                    for(String seq: this.getAdapterSequences()){
                        cmd.add("-A " + seq);
                    }
                }
                
                if(this.getIncludeRefSeq())
                    cmd.add("-R");
                if(!this.getSkipSAMHeader())
                    cmd.add("-H");
                if(this.getReportUnmappedReads())
                    cmd.add("-u");
                
                cmd.add("-B " + this.getStartAtThisRead());
                cmd.add("-H " + this.getEndtAtThisRead());
                cmd.add("-D " + this.getDigestionSite());
                cmd.add("-S " + this.getSeedSize());
                
                /*
                    -n 0: only map to 2 forward strands, i.e. BSW(++) and BSC(-+) ("Lister protocol")
                    for PE sequencing, map read#1 to ++ and -+, read#2 to +- and --.
                    -n 1: map SE or PE reads to all 4 strands, i.e. ++, +-, -+, -- ("Cokus protocol")
                    default: -n 0. Most bisulfite sequencing data is generated only from forward strands.
                */
                cmd.add("-n " + this.getMappingStrand());
                
                /*
                  set the alignment information for the additional nucleotide transition. <str> is in the form of two different nucleotides,
                  the first one in the reads could be mapped to the second one in the reference sequences.
                  default: -M TC, corresponds to C=>U(T) transition in bisulfite conversion.
                  example: -M GA could be used to detect to A=>I(G) transition in RNA editing.
                */
                cmd.add("-M " + this.getTransitionMap());
                
                /*
                verbose level:  0=no message displayed (quiet mode); 
                                1=major message (default); 
                                2=detailed message                
                */
                cmd.add("-V " + this.getMessageLevel());
                
                /*
                  still not sure how to handle these two parameters
                    -w  <int>   max number of equal best hits to count, smaller will be faster, default=MAXHITS in makefile
                    -z  <int>   base quality, default=33 [Illumina is using 64, Sanger Institute is using 33]
                */                
                                
            }
            catch(Exception ex){
                getLogger().info("error executing bsmap command:\n" + ex.toString());
                getLogger().error("error executing bsmap command:\n" + ex.toString());
                throw new IOException("error executing bsmap command");
            }
        }
        
        getLogger().info(STEP_ID_STRING + ": completed");
    }
    
    
    
            
    /**
     * this should be called prior to executing the step.
     * check unzip software exists and input files are available
     * 
     * @throws IOException
     */
    @Override
    public void verifyInputData() throws IOException{
        
        getLogger().info("verify input data");        
        this.setPaths();
        
        if(new File(this.getPathToSoftware()).exists() == false){
            getLogger().info("BSMap software not found at location < " + this.getPathToSoftware() +">");
            getLogger().error("BSMap software not found at location < " + this.getPathToSoftware() +">");
            throw new IOException("unzip software not found at location < " + this.getPathToSoftware() +">");
        }
                
                
                            
        // check the data files
        String fastqFile1in = "";
        String fastqFile1out = "";
        Iterator itSD = this.stepInputData.getSampleData().iterator();
        for(SampleDataEntry sampleData: this.stepInputData.getSampleData()){
            
            
            //Fastq 1
            if (sampleData.getFastqFile1()==null) {
                getLogger().error("no Fastq1 file specified");
                throw new IOException("no Fastq1 file specified");
            }
            
            fastqFile1out = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1());
            fastqFile1in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION));
            String fastqFile1 = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION);
            if (new File(fastqFile1in).exists()==false && new File(fastqFile1out).exists()==false){
                getLogger().error(STEP_ID_STRING + ": fastq 1 Files <" + fastqFile1in + "> & <" + fastqFile1out + "> do not exist");
                throw new IOException(STEP_ID_STRING + ": fastq 2 Files <" + fastqFile1in + "> & <" + fastqFile1out + "> do not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                getLogger().info(STEP_ID_STRING + ": incorrect file extension for input file <" 
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
                getLogger().error(STEP_ID_STRING + ": fastq 2 Files <" + fastqFile2in + "> & <" + fastqFile2out + "> do not exist");
                throw new IOException(STEP_ID_STRING + ": fastq 2 Files <" + fastqFile2in + "> & <" + fastqFile2out + "> do not exist");
            }
            if (fastqFile2in.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                getLogger().error(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                  + fastqFile2in + ">. should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                  + fastqFile2in + ">. \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
                        
            
        }

    }
    
    
    
    @Override
    public void verifyOutputData(){
        getLogger().info("no output verification required");
    }
    
    
    
    /**
     * generate sample configuration data so the user can see what can be 
     * specified 
     * 
     * @return 
     */
    @Override
    public HashMap generateExampleConfigurationData() {

        getLogger().info(STEP_ID_STRING + ": generate example configuration data");
        
        HashMap<String, Object> configData = new HashMap();
        
        configData.put(ID_SOFTWARE, "/usr/local/bin/bsmap");
        configData.put(ID_REF_GENOME, "hsa");
        configData.put(ID_SEED_SIZE, 16);
        configData.put(ID_MISMATCHES, 0.08);
        configData.put(ID_GAPSIZE, 0);
        configData.put(ID_THREEMAPPING, false);
        configData.put(ID_THREADS, 4);
        configData.put(ID_KMERCUTOFF, 1e-06);
        configData.put(ID_REPORT_REPEATS, 1);
        configData.put(ID_QUAL_THRESHOLD_VALUE, 0);
        configData.put(ID_LOW_QUAL_VALUE, 5);
        configData.put(ID_MAX_INSERT_SIZE, 500);
        configData.put(ID_MIN_INSERT_SIZE, 28);
        configData.put(ID_MAP_FIRST_N_NUCS, 0);
        configData.put(ID_GENOME_INDEX_INTERVAL, 4);
        configData.put(ID_TRIM_ADAPTERSEQ, true);
        configData.put(ID_ADAPTERSEQ, new ArrayList<>(Arrays.asList(
                "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG", 
                "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")));        
        configData.put(ID_ADAPTERSEQ, "");
        configData.put(ID_INCLUDE_REF_SEQ, false);
        configData.put(ID_SKIP_SAM_HEADER, false);
        configData.put(ID_REPORT_UNMAPPED_READS, false);
        configData.put(ID_START_AT_READ_N, 1);
        configData.put(ID_END_AT_READ_N, -1);
        configData.put(ID_DIGESTIONSITE, "C-CGG");
        configData.put(ID_RANDOMSEED, 0);
        configData.put(ID_MAPPINGSTRAND, 0);
        configData.put(ID_TRANSITIONMAP, "TC");
        configData.put(ID_MESSAGELEVEL, 1);

        return configData;
        
    }

    
    
    
    
    
    /**
     * @return the unzipSoftware
     */
    public String getPathToSoftware() {
        return pathToSoftware;
    }

    /**
     * @param unzipSoftware the unzipSoftware to set
     */
    public void setPathToSoftware(String unzipSoftware) {
        this.pathToSoftware = unzipSoftware;
    }

    /**
     * @return the pathToBSMap
     */
    public String getPathToBSMap() {
        return pathToBSMap;
    }

    /**
     * @param pathToBSMap the pathToBSMap to set
     */
    public void setPathToBSMap(String pathToBSMap) {
        this.pathToBSMap = pathToBSMap;
    }

    /**
     * @return the seedSize
     */
    public int getSeedSize() {
        return seedSize;
    }

    /**
     * @param seedSize the seedSize to set
     */
    public void setSeedSize(int seedSize) {
        this.seedSize = seedSize;
    }

    /**
     * @return the AlignMode
     */
    public String getAlignMode() {
        return AlignMode;
    }

    /**
     * @param AlignMode the AlignMode to set
     */
    public void setAlignMode(String AlignMode) {
        this.AlignMode = AlignMode;
    }

    /**
     * @return the noOfMismatches
     */
    public double getNoOfMismatches() {
        return noOfMismatches;
    }

    /**
     * @param noOfMismatches the noOfMismatches to set
     */
    public void setNoOfMismatches(double noOfMismatches) {
        this.noOfMismatches = noOfMismatches;
    }

    /**
     * @return the threeNucMapping
     */
    public Boolean getThreeNucMapping() {
        return threeNucMapping;
    }

    /**
     * @param threeNucMapping the threeNucMapping to set
     */
    public void setThreeNucMapping(Boolean threeNucMapping) {
        this.threeNucMapping = threeNucMapping;
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
     * @return the rootDataFolder
     */
    public String getRootDataFolder() {
        return rootDataFolder;
    }

    /**
     * @param rootDataFolder the rootDataFolder to set
     */
    public void setRootDataFolder(String rootDataFolder) {
        this.rootDataFolder = rootDataFolder;
    }

    /**
     * @return the ReferenceGenome
     */
    public String getReferenceGenome() {
        return ReferenceGenome;
    }

    /**
     * @param ReferenceGenome the ReferenceGenome to set
     */
    public void setReferenceGenome(String ReferenceGenome) {
        this.ReferenceGenome = ReferenceGenome;
    }

    /**
     * @return the kmerCutoffRatio
     */
    public double getKmerCutoffRatio() {
        return kmerCutoffRatio;
    }

    /**
     * @param kmerCutoffRatio the kmerCutoffRatio to set
     */
    public void setKmerCutoffRatio(double kmerCutoffRatio) {
        this.kmerCutoffRatio = kmerCutoffRatio;
    }

    /**
     * @return the reportRepeatsValue
     */
    public int getReportRepeatsValue() {
        return reportRepeatsValue;
    }

    /**
     * @param reportRepeatsValue the reportRepeatsValue to set
     */
    public void setReportRepeatsValue(int reportRepeatsValue) {
        this.reportRepeatsValue = reportRepeatsValue;
    }

    /**
     * @return the qualityThresholdTrimValue
     */
    public int getQualityThresholdTrimValue() {
        return qualityThresholdTrimValue;
    }

    /**
     * @param qualityThresholdTrimValue the qualityThresholdTrimValue to set
     */
    public void setQualityThresholdTrimValue(int qualityThresholdTrimValue) {
        this.qualityThresholdTrimValue = qualityThresholdTrimValue;
    }

    /**
     * @return the lowQualityFilterValue
     */
    public int getLowQualityFilterValue() {
        return lowQualityFilterValue;
    }

    /**
     * @param lowQualityFilterValue the lowQualityFilterValue to set
     */
    public void setLowQualityFilterValue(int lowQualityFilterValue) {
        this.lowQualityFilterValue = lowQualityFilterValue;
    }

    /**
     * @return the maxInsertSize
     */
    public int getMaxInsertSize() {
        return maxInsertSize;
    }

    /**
     * @param maxInsertSize the maxInsertSize to set
     */
    public void setMaxInsertSize(int maxInsertSize) {
        this.maxInsertSize = maxInsertSize;
    }

    /**
     * @return the minInsertSize
     */
    public int getMinInsertSize() {
        return minInsertSize;
    }

    /**
     * @param minInsertSize the minInsertSize to set
     */
    public void setMinInsertSize(int minInsertSize) {
        this.minInsertSize = minInsertSize;
    }

    /**
     * @return the mapFirstNucleotides
     */
    public int getMapFirstNucleotides() {
        return mapFirstNucleotides;
    }

    /**
     * @param mapFirstNucleotides the mapFirstNucleotides to set
     */
    public void setMapFirstNucleotides(int mapFirstNucleotides) {
        this.mapFirstNucleotides = mapFirstNucleotides;
    }

    /**
     * @return the genomeIndexInterval
     */
    public int getGenomeIndexInterval() {
        return genomeIndexInterval;
    }

    /**
     * @param genomeIndexInterval the genomeIndexInterval to set
     */
    public void setGenomeIndexInterval(int genomeIndexInterval) {
        this.genomeIndexInterval = genomeIndexInterval;
    }

    /**
     * @return the includeRefSeq
     */
    public Boolean getIncludeRefSeq() {
        return includeRefSeq;
    }

    /**
     * @param includeRefSeq the includeRefSeq to set
     */
    public void setIncludeRefSeq(Boolean includeRefSeq) {
        this.includeRefSeq = includeRefSeq;
    }

    /**
     * @return the skipSAMHeader
     */
    public Boolean getSkipSAMHeader() {
        return skipSAMHeader;
    }

    /**
     * @param skipSAMHeader the skipSAMHeader to set
     */
    public void setSkipSAMHeader(Boolean skipSAMHeader) {
        this.skipSAMHeader = skipSAMHeader;
    }

    /**
     * @return the reportUnmappedReads
     */
    public Boolean getReportUnmappedReads() {
        return reportUnmappedReads;
    }

    /**
     * @param reportUnmappedReads the reportUnmappedReads to set
     */
    public void setReportUnmappedReads(Boolean reportUnmappedReads) {
        this.reportUnmappedReads = reportUnmappedReads;
    }

    /**
     * @return the startAtThisRead
     */
    public int getStartAtThisRead() {
        return startAtThisRead;
    }

    /**
     * @param startAtThisRead the startAtThisRead to set
     */
    public void setStartAtThisRead(int startAtThisRead) {
        this.startAtThisRead = startAtThisRead;
    }

    /**
     * @return the endtAtThisRead
     */
    public int getEndtAtThisRead() {
        return endtAtThisRead;
    }

    /**
     * @param endtAtThisRead the endtAtThisRead to set
     */
    public void setEndtAtThisRead(int endtAtThisRead) {
        this.endtAtThisRead = endtAtThisRead;
    }

    /**
     * @return the digestionSite
     */
    public String getDigestionSite() {
        return digestionSite;
    }

    /**
     * @param digestionSite the digestionSite to set
     */
    public void setDigestionSite(String digestionSite) {
        this.digestionSite = digestionSite;
    }

    /**
     * @return the randomSeed
     */
    public int getRandomSeed() {
        return randomSeed;
    }

    /**
     * @param randomSeed the randomSeed to set
     */
    public void setRandomSeed(int randomSeed) {
        this.randomSeed = randomSeed;
    }

    /**
     * @return the mappingStrand
     */
    public int getMappingStrand() {
        return mappingStrand;
    }

    /**
     * @param mappingStrand the mappingStrand to set
     */
    public void setMappingStrand(int mappingStrand) {
        this.mappingStrand = mappingStrand;
    }

    /**
     * @return the transitionMap
     */
    public String getTransitionMap() {
        return transitionMap;
    }

    /**
     * @param transitionMap the transitionMap to set
     */
    public void setTransitionMap(String transitionMap) {
        this.transitionMap = transitionMap;
    }

    /**
     * @return the gapSize
     */
    public int getGapSize() {
        return gapSize;
    }

    /**
     * @param gapSize the gapSize to set
     */
    public void setGapSize(int gapSize) {
        this.gapSize = gapSize;
    }

    /**
     * @return the messageLevel
     */
    public int getMessageLevel() {
        return messageLevel;
    }

    /**
     * @param messageLevel the messageLevel to set
     */
    public void setMessageLevel(int messageLevel) {
        this.messageLevel = messageLevel;
    }

    /**
     * @return the trimAdapterSequence
     */
    public Boolean getTrimAdapterSequence() {
        return trimAdapterSequence;
    }

    /**
     * @param trimAdapterSequence the trimAdapterSequence to set
     */
    public void setTrimAdapterSequence(Boolean trimAdapterSequence) {
        this.trimAdapterSequence = trimAdapterSequence;
    }

    /**
     * @return the adapterSequences
     */
    public ArrayList<String> getAdapterSequences() {
        return adapterSequences;
    }

    /**
     * @param adapterSequences the adapterSequences to set
     */
    public void setAdapterSequences(ArrayList<String> adapterSequences) {
        this.adapterSequences = adapterSequences;
    }
    
    
    
    
    
}
