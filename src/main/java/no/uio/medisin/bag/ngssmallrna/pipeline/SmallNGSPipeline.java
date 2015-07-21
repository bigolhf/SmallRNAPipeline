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
    
    private String                      trimAdapterFile = "";
    private int                         trimNoOfMismatches = 2;
    private int                         trimMinAlignScore = 7;
    private int                         trimNoOfThreads = 4;
    private int                         trimMinAvgReadQuality = 30;
    
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
                    
                    StepInputData sid = new StepInputData(trimAdapterParams, this.getPipelineData().getProjectID(), this.getPipelineData().getProjectRoot(), this.getSampleData());
                    TrimAdaptersStep ngsStep = new TrimAdaptersStep(sid);
                    this.getNGSSteps().add(ngsStep);
                    ngsStep.verifyInputData();
                    ngsStep.execute();
                    
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

        HashMap softwareOptions = (HashMap) pipelineConfiguration.get("software");
        this.setSoftwareRootFolder((String) softwareOptions.get("root_folder"));
        this.setAdapterTrimmingSoftware((String) softwareOptions.get("adapter_trimming"));
        
        HashMap trimAdapterOptions = (HashMap) pipelineConfiguration.get("adapter_trimming");
        this.setTrimAdapterFile((String) trimAdapterOptions.get("adapter_file"));
        this.setTrimMinAlignScore((int) trimAdapterOptions.get("no_of_mismatches"));
        this.setTrimNoOfMismatches((int) trimAdapterOptions.get("min_align_score"));
        this.setTrimNoOfThreads((int) trimAdapterOptions.get("no_of_threads"));
        this.setTrimMinAvgReadQuality((int) trimAdapterOptions.get("min_avg_read_qual"));
        
                        
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
    
    
    
    
}
