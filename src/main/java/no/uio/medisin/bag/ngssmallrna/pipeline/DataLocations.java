/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

/**
 * contains information about root paths to the different data types
 * 
 * @author simonray
 */
public class DataLocations {
    
    public final static     String ID_REL_BOWTIE_PATH       = "Sequence/BowtieIndex/genome";
    public final static     String ID_REL_ABUN_DATA_PATH    = "/Sequence/AbundantSequences/abundant";
    public final static     String ID_REL_WHOLE_GENSEQ_PATH = "Sequence/WholeGenomeFasta";
    public final static     String ID_GENE_ANNOTATION       = "/Annotation/Genes";
    
    public final static     String ID_GENOME_FOLDER = "genome_root_folder:";
    public final static     String ID_MIRBASE_FOLDER = "mirbase_folder:";
    
    
    private       String genomeRootFolder;
    private       String mirbaseFolder;

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
    
}
