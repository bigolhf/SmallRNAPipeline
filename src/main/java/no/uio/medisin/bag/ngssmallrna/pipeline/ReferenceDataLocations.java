/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import static no.uio.medisin.bag.ngssmallrna.pipeline.MiRNAFeature.logger;

/**
 * contains information about root paths to the different data types
 * 
 * @author simonray
 */
public class ReferenceDataLocations {
    
    public final static     String ID_CONFIG_ID             = "ReferenceData";
    public final static     String ID_GENOME_FOLDER         = "genomeRootFolder";
    public final static     String ID_MIRBASE_FOLDER        = "mirbaseFolder";
    public final static     String ID_TARGETSCAN_FOLDER     = "targetscanFolder";
    
    public final static     String ID_REL_BOWTIE_PATH       = "Sequence/BowtieIndex/genome";
    public final static     String ID_REL_ABUN_DATA_PATH    = "/Sequence/AbundantSequences/abundant";
    public final static     String ID_REL_WHOLE_GENSEQ_PATH = "Sequence/WholeGenomeFasta";
    public final static     String ID_GENE_ANNOTATION       = "/Annotation/Genes";
    
    
    
    private                 String genomeRootFolder;
    private                 String mirbaseFolder;
    private                 String targetbaseFolder;

    
    public ReferenceDataLocations(){
    
    }
    
    
    public ReferenceDataLocations(HashMap options){
        
        genomeRootFolder    = (String) options.get(ReferenceDataLocations.ID_GENOME_FOLDER);
        mirbaseFolder       = (String)options.get(ReferenceDataLocations.ID_MIRBASE_FOLDER);
        
    }
    /**
     * in this method we are simply checking that the configuration file 
     * has all the entries we need. We dont check if the values are acceptable
     * that is the role of the NGSStep.
     * 
     * @param configData
     * @throws Exception 
     */
    
    public void parseConfigurationData(HashMap configData) throws Exception{
        logger.info(this.getClass() + ": verify configuration data");
        
        if(configData.get(ID_GENOME_FOLDER)==null) {
            throw new NullPointerException("<" + ID_GENOME_FOLDER + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIRBASE_FOLDER)==null) {
            throw new NullPointerException("<" + ID_MIRBASE_FOLDER + "> : Missing Definition in Configuration File");
        }

        logger.info("passed");
    }
    
    
    
    /**
     * 
     * @throws IOException 
     */
    public void verifyReferenceData() throws IOException{
        if(new File(ID_GENOME_FOLDER).exists()==false){
            throw new IOException("root genome folder <" + ID_GENOME_FOLDER + "> not found");
        }
        if(new File(ID_MIRBASE_FOLDER).exists()==false){
            throw new IOException("root mirbase folder <" + ID_MIRBASE_FOLDER + "> not found");
        }
    }
    

    
    /**
     * generate sample configuration data
     * 
     * @return 
     */
    public HashMap generateExampleConfigurationData() {

        HashMap configData = new HashMap();
        configData.put(ID_GENOME_FOLDER, "/data/genomes");
        configData.put(ID_MIRBASE_FOLDER, "/data/mirbase");
    
        return configData;
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

    /**
     * @return the mirbaseFolder
     */
    public String getMirbaseFolder() {
        return mirbaseFolder;
    }

    /**
     * @param mirbaseFolder the mirbaseFolder to set
     */
    public void setMirbaseFolder(String mirbaseFolder) {
        this.mirbaseFolder = mirbaseFolder;
    }

    /**
     * @return the targetbaseFolder
     */
    public String getTargetbaseFolder() {
        return targetbaseFolder;
    }

    /**
     * @param targetbaseFolder the targetbaseFolder to set
     */
    public void setTargetbaseFolder(String targetbaseFolder) {
        this.targetbaseFolder = targetbaseFolder;
    }
    
}
