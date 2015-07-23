/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import static no.uio.medisin.bag.ngssmallrna.pipeline.SmallNGSCmd.logger;
import no.uio.medisin.bag.ngssmallrna.steps.BowtieMapReadsStep;
import no.uio.medisin.bag.ngssmallrna.steps.CollapseReadsStep;
import no.uio.medisin.bag.ngssmallrna.steps.NGSStep;
import no.uio.medisin.bag.ngssmallrna.steps.NGSRunStepData;
import no.uio.medisin.bag.ngssmallrna.steps.StepInputData;
import no.uio.medisin.bag.ngssmallrna.steps.TrimAdaptersStep;
import org.yaml.snakeyaml.Yaml;
import org.codehaus.jackson.JsonGenerationException;
import org.codehaus.jackson.annotate.JsonAutoDetect;
import org.codehaus.jackson.annotate.JsonMethod;
import org.codehaus.jackson.map.JsonMappingException;
import org.codehaus.jackson.map.ObjectMapper;

/**
 *
 * @author sr
 */
public class SmallNGSPipeline {
    
    static  Yaml                        yaml = new Yaml();
    static  String                      FileSeparator = System.getProperty("file.separator");
    
    private String                      configurationFile = "";
    private String                      pipelineFile = "";
    private String                      dataFile = "";
    
    private String                      softwareRootFolder = "";
    private String                      adapterTrimmingSoftware = "";
    private String                      fastq2fastaSoftware = "";
    private String                      collapseFastaSoftware = "";
    
    private String                      genomeRootFolder = "";
    
    private String                      trimAdapterFile = "";
    private int                         trimNoOfMismatches = 2;
    private int                         trimMinAlignScore = 7;
    private int                         trimNoOfThreads = 4;
    private int                         trimMinAvgReadQuality = 30;
    
    private String                      bowtieMappingAlignMode = "";
    private int                         bowtieMappingNoOfMismatches = 2;
    private int                         bowtieMappingNoOfThreads = 4;
    private String                      bowtieMappingReferenceGenome = "";
    
    private PipelineData                pipelineData = new PipelineData();
    private ArrayList<SampleDataEntry>  SampleData = new ArrayList<>();
    private ArrayList<NGSStep>          NGSSteps = new ArrayList<>();

    
    
    /**
     * here we load the information from the three specified files
     * (run, pipeline and dataset) 
     * 
     * @throws IOException 
     */
    public void prepare_pipeline() throws IOException {      
        
        
        
        if(new File(this.getConfigurationFile()).exists()== false)
        {
            throw new IOException("configuration file <" + this.getConfigurationFile() + "> does not exist");
        }
        
        if(new File(this.getPipelineFile()).exists()== false)
        {
            throw new IOException("pipeline file <" + this.getPipelineFile() + "> does not exist");
        }
        
        if(new File(this.getDataFile()).exists()== false)
        {
            throw new IOException("data file <" + this.getDataFile() + "> does not exist");
        }
        
        this.readConfigurationFile();
        this.readPipelineFile();
        this.readDataFile();
        
        
    }
    
    
    /**
     * 
     * build each step and execute
     * 
     */
    public void buildPipeline(){
        
        for (NGSRunStepData stepData: this.getPipelineData().getStepsData()){
            
            switch (stepData.getStepType()){
                case "TrimAdapters":
                    
                    HashMap trimAdapterParams = new HashMap();
                    trimAdapterParams.put("trimAdapterSoftware",    this.getSoftwareRootFolder() + FileSeparator + this.getAdapterTrimmingSoftware());
                    trimAdapterParams.put("trimAdapterFile",        this.getTrimAdapterFile());
                    trimAdapterParams.put("trimNoOfMismatches",     this.getTrimNoOfMismatches());
                    trimAdapterParams.put("trimMinAlignScore",      this.getTrimMinAlignScore());
                    trimAdapterParams.put("trimNoOfThreads",        this.getTrimNoOfThreads());
                    trimAdapterParams.put("trimMinAvgReadQuality",  this.getTrimMinAvgReadQuality());
                    
                    StepInputData sidTrim = new StepInputData(trimAdapterParams, this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), this.getSampleData());
                    TrimAdaptersStep ngsTrimStep = new TrimAdaptersStep(sidTrim);
                    this.getNGSSteps().add(ngsTrimStep);
                    ngsTrimStep.verifyInputData();
                    ngsTrimStep.execute();
                    
                    break;
                    
                    
                    
                case "CollapseReads":
                    
                    HashMap collapseReadsParams = new HashMap();
                    collapseReadsParams.put("fastqTofasta", this.getFastq2fastaSoftware());
                    collapseReadsParams.put("collapseFasta", this.getCollapseFastaSoftware());
                    
                    StepInputData sidCollapse = new StepInputData(collapseReadsParams, this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), this.getSampleData());
                    CollapseReadsStep ngsCollapseStep = new CollapseReadsStep(sidCollapse);
                    ngsCollapseStep.verifyInputData();
                    ngsCollapseStep.execute();
                    
                    break;
                    
                case "BowtieMapReads":
                    
                    HashMap bowtieMapReadsParams = new HashMap();
                    bowtieMapReadsParams.put("bowtieMapGenomeRootFolder", this.getGenomeRootFolder());
                    bowtieMapReadsParams.put("bowtieReferenceGenome", this.getBowtieMappingReferenceGenome());
                    bowtieMapReadsParams.put("bowtieMapAlignMode", this.getBowtieMappingAlignMode());
                    bowtieMapReadsParams.put("bowtieMapNoOfMismatches", this.getBowtieMappingNoOfMismatches()); //getBowtieMappingNoOfMismatches
                    bowtieMapReadsParams.put("bowtieNoOfThreads", this.getBowtieMappingNoOfThreads());
                    
                    StepInputData sidMap = new StepInputData(bowtieMapReadsParams, this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), this.getSampleData());
                    BowtieMapReadsStep ngsBowtieMapReads = new BowtieMapReadsStep(sidMap);
                    ngsBowtieMapReads.verifyInputData();
                    ngsBowtieMapReads.execute();
                    
                    break;
                    
            }
            
        }
        
    }

    

    
    public void runPipeline(){
        
        
        for (NGSStep ngsStep: this.getNGSSteps()){
            
        }
        
    }
    
    
    /**
     * read file containing information about the data to be analyzed
     * 
     * @throws IOException 
     * 
     */
    public void readDataFile() throws IOException{
        
        
        logger.info("read data file");
        String line = "";
        BufferedReader bwData = new BufferedReader(new FileReader(new File(this.getDataFile())));
            bwData.readLine(); // skip header line
            while((line=bwData.readLine())!= null){
                
                String tokens[] = line.split("\t");
                String file = tokens[0];
                String source = tokens[1];
                String condition = tokens[2];
                String time = tokens[3];
                
                getSampleData().add(new SampleDataEntry(file, source, condition, time));
                
            }
        bwData.close();
        logger.info("read " + getSampleData().size() + " entries");
        
        
    }
    
    
    
    
    /**
     * read configuration settings for the pipeline
     * 
     * @throws IOException 
     */
    public void readConfigurationFile() throws IOException{
        
        logger.info("read pipeline configuration file");
        HashMap pipelineConfiguration = (HashMap) yaml.load(new FileInputStream(new File(getConfigurationFile())));
        
        HashMap dataOptions = (HashMap) pipelineConfiguration.get("data");
        this.setGenomeRootFolder((String) dataOptions.get("genome_root_folder"));

        HashMap softwareOptions = (HashMap) pipelineConfiguration.get("software");
        this.setSoftwareRootFolder((String) softwareOptions.get("root_folder"));
        this.setAdapterTrimmingSoftware((String) softwareOptions.get("adapter_trimming"));
        this.setFastq2fastaSoftware((String) softwareOptions.get("fastq_to_fasta"));
        this.setCollapseFastaSoftware((String) softwareOptions.get("fastx_collapser"));
        
        HashMap trimAdapterOptions = (HashMap) pipelineConfiguration.get("adapter_trimming");
        this.setTrimAdapterFile((String) trimAdapterOptions.get("adapter_file"));
        this.setTrimMinAlignScore((int) trimAdapterOptions.get("no_of_mismatches"));
        this.setTrimNoOfMismatches((int) trimAdapterOptions.get("min_align_score"));
        this.setTrimNoOfThreads((int) trimAdapterOptions.get("no_of_threads"));
        this.setTrimMinAvgReadQuality((int) trimAdapterOptions.get("min_avg_read_qual"));
        
        HashMap mapReadsOptions = (HashMap) pipelineConfiguration.get("bowtie_mapping");
        this.setBowtieMappingAlignMode((String) mapReadsOptions.get("alignment_mode"));
        this.setBowtieMappingNoOfMismatches((int) mapReadsOptions.get("no_of_mismatches"));
        this.setBowtieMappingNoOfThreads((int) mapReadsOptions.get("no_of_threads"));
        this.setBowtieMappingReferenceGenome((String) mapReadsOptions.get("host"));
        
        logger.info("done\n");
        
    }
    
    
    
    
    /**
     * 
     */
    public void readPipelineFile(){
        
	ObjectMapper mapper = new ObjectMapper();
        mapper.setVisibility(JsonMethod.FIELD, JsonAutoDetect.Visibility.ANY);
        
	try { 
            
            
                logger.info("read pipeline file <" + this.getPipelineFile() + ">");
                pipelineData = mapper.readValue(new File(this.getPipelineFile()), PipelineData.class); 
                logger.info("done");
                logger.info(getPipelineData().toString());
                
 
	} catch (JsonGenerationException e) {
 
		e.printStackTrace();
 
	} catch (JsonMappingException e) {
 
		e.printStackTrace();
 
	} catch (IOException e) {
 
		e.printStackTrace();
 
	}
        
        
    }
    
    
    
    
    

    /**
     * @return the configurationFile
     */
    public String getConfigurationFile() {
        return configurationFile;
    }

    /**
     * @param configurationFile the configurationFile to set
     */
    public void setConfigurationFile(String configurationFile) {
        this.configurationFile = configurationFile;
    }

    /**
     * @return the pipelineFile
     */
    public String getPipelineFile() {
        return pipelineFile;
    }

    /**
     * @param pipelineFile the pipelineFile to set
     */
    public void setPipelineFile(String pipelineFile) {
        this.pipelineFile = pipelineFile;
    }

    /**
     * @return the dataFile
     */
    public String getDataFile() {
        return dataFile;
    }

    /**
     * @param dataFile the dataFile to set
     */
    public void setDataFile(String dataFile) {
        this.dataFile = dataFile;
    }

    /**
     * @return the pipelineData
     */
    public PipelineData getPipelineData() {
        return pipelineData;
    }

    /**
     * @return the softwareRootFolder
     */
    public String getSoftwareRootFolder() {
        return softwareRootFolder;
    }

    /**
     * @param softwareRootFolder the softwareRootFolder to set
     */
    public void setSoftwareRootFolder(String softwareRootFolder) {
        this.softwareRootFolder = softwareRootFolder;
    }

    /**
     * @return the adapterTrimmingSoftware
     */
    public String getAdapterTrimmingSoftware() {
        return adapterTrimmingSoftware;
    }

    /**
     * @param adapterTrimmingSoftware the adapterTrimmingSoftware to set
     */
    public void setAdapterTrimmingSoftware(String adapterTrimmingSoftware) {
        this.adapterTrimmingSoftware = adapterTrimmingSoftware;
    }

    /**
     * @return the NGSSteps
     */
    public ArrayList<NGSStep> getNGSSteps() {
        return NGSSteps;
    }

    /**
     * @return the trimAdapterFile
     */
    public String getTrimAdapterFile() {
        return trimAdapterFile;
    }

    /**
     * @param trimAdapterFile the trimAdapterFile to set
     */
    public void setTrimAdapterFile(String trimAdapterFile) {
        this.trimAdapterFile = trimAdapterFile;
    }

    /**
     * @return the trimNoOfMismatches
     */
    public int getTrimNoOfMismatches() {
        return trimNoOfMismatches;
    }

    /**
     * @param trimNoOfMismatches the trimNoOfMismatches to set
     */
    public void setTrimNoOfMismatches(int trimNoOfMismatches) {
        this.trimNoOfMismatches = trimNoOfMismatches;
    }

    /**
     * @return the trimMinAlignScore
     */
    public int getTrimMinAlignScore() {
        return trimMinAlignScore;
    }

    /**
     * @param trimMinAlignScore the trimMinAlignScore to set
     */
    public void setTrimMinAlignScore(int trimMinAlignScore) {
        this.trimMinAlignScore = trimMinAlignScore;
    }

    /**
     * @return the trimNoOfThreads
     */
    public int getTrimNoOfThreads() {
        return trimNoOfThreads;
    }

    /**
     * @param trimNoOfThreads the trimNoOfThreads to set
     */
    public void setTrimNoOfThreads(int trimNoOfThreads) {
        this.trimNoOfThreads = trimNoOfThreads;
    }

    /**
     * @return the SampleData
     */
    public ArrayList<SampleDataEntry> getSampleData() {
        return SampleData;
    }

    /**
     * @return the trimMinReadQuality
     */
    public int getTrimMinAvgReadQuality() {
        return trimMinAvgReadQuality;
    }

    /**
     * @param trimMinReadQuality the trimMinReadQuality to set
     */
    public void setTrimMinAvgReadQuality(int trimMinReadQuality) {
        this.trimMinAvgReadQuality = trimMinReadQuality;
    }

    /**
     * @return the fastq2fastaSoftware
     */
    public String getFastq2fastaSoftware() {
        return fastq2fastaSoftware;
    }

    /**
     * @param fastq2fastaSoftware the fastq2fastaSoftware to set
     */
    public void setFastq2fastaSoftware(String fastq2fastaSoftware) {
        this.fastq2fastaSoftware = fastq2fastaSoftware;
    }

    /**
     * @return the collapseFastaSoftware
     */
    public String getCollapseFastaSoftware() {
        return collapseFastaSoftware;
    }

    /**
     * @param collapseFastaSoftware the collapseFastaSoftware to set
     */
    public void setCollapseFastaSoftware(String collapseFastaSoftware) {
        this.collapseFastaSoftware = collapseFastaSoftware;
    }

    /**
     * @return the bowtieMappingAlignMode
     */
    public String getBowtieMappingAlignMode() {
        return bowtieMappingAlignMode;
    }

    /**
     * @param bowtieMappingAlignMode the bowtieMappingAlignMode to set
     */
    public void setBowtieMappingAlignMode(String bowtieMappingAlignMode) {
        this.bowtieMappingAlignMode = bowtieMappingAlignMode;
    }

    /**
     * @return the bowtieMappingNoOfMismatches
     */
    public int getBowtieMappingNoOfMismatches() {
        return bowtieMappingNoOfMismatches;
    }

    /**
     * @param bowtieMappingNoOfMismatches the bowtieMappingNoOfMismatches to set
     */
    public void setBowtieMappingNoOfMismatches(int bowtieMappingNoOfMismatches) {
        this.bowtieMappingNoOfMismatches = bowtieMappingNoOfMismatches;
    }

    /**
     * @return the bowtieMappingNoOfThreads
     */
    public int getBowtieMappingNoOfThreads() {
        return bowtieMappingNoOfThreads;
    }

    /**
     * @param bowtieMappingNoOfThreads the bowtieMappingNoOfThreads to set
     */
    public void setBowtieMappingNoOfThreads(int bowtieMappingNoOfThreads) {
        this.bowtieMappingNoOfThreads = bowtieMappingNoOfThreads;
    }

    /**
     * @return the bowtieMappingReferenceGenome
     */
    public String getBowtieMappingReferenceGenome() {
        return bowtieMappingReferenceGenome;
    }

    /**
     * @param bowtieMappingReferenceGenome the bowtieMappingReferenceGenome to set
     */
    public void setBowtieMappingReferenceGenome(String bowtieMappingReferenceGenome) {
        this.bowtieMappingReferenceGenome = bowtieMappingReferenceGenome;
    }

    /**
     * @return the genomeRootFolder
     */
    public String getGenomeRootFolder() {
        return genomeRootFolder;
    }

    /**
     * @param genomeRootFolder the genomeRootFolder to set
     */
    public void setGenomeRootFolder(String genomeRootFolder) {
        this.genomeRootFolder = genomeRootFolder;
    }
    
    
    
    
}
