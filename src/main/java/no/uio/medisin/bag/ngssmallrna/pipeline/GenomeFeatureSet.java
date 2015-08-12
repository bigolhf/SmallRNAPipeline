/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

/**
 * handles access and searching a set
 * @author sr
 */
public class GenomeFeatureSet {
    
    static Logger logger = LogManager.getLogger();
    
    private List<GenomeFeature>          genomeFeatureList;

    
    
    public GenomeFeatureSet(){
        genomeFeatureList = new ArrayList<>();
    }
    
    
    
    /**
     * Does the read sufficiently overlap a defined Feature entry?
     * 
     * @param start
     * @param stop
     * @param chr
     * @param strand
     * @param bleed         int : specifies how much a read can 'miss' an entry
     *                            and still be counted
     * @param featureTypes  String : search only features of this type
     * 
     * @return MiRNAFeature
     */
    public GenomeFeature doesReadOverlapFeature(int start, int stop, String chr, String strand, int bleed, String featureTypes[]){
        
        for(GenomeFeature feature: this.genomeFeatureList){
            
            if (featureTypes[0].equals("all") || Arrays.asList(featureTypes).contains(feature.getFeatureName())){
                
                if (feature.chromosomeMatch(chr)){
                    
                    if(strand.equals(feature.getStrand())){
                        if((start - bleed) <= feature.getStartPos()){

                            if( feature.getEndPos() <= (stop + bleed)){
                                logger.info("\t" + start + "\t" + stop + "\t" + feature.getStartPos());
                                return feature;
                            }

                        }

                    }                
                }
                
            }
            
        }
        
        return null;
        
    }

    
    
    
    
    /**
     * add feature to list
     * 
     * @param gffString 
     */
    public void addFeature(String gffString){
        if (gffString.contains("exon"))
            genomeFeatureList.add(new GenomeFeature(gffString));
    }
    
    
    
    /**
     * 1. 
     * load miRNA specs (name and Chromosome position) from GFF file
     * downloaded from miRBase
     * Because different releases of miRBase use different releases of
     * reference genome, we have to track both miRBase and genome reference IDs
     * 2.
     * Load Sequence from Mature.fa file
     * 
     * @param host              String : 3 char abbreviation for host
     * @param genomeGFFFile    String : absolute path to file
     * @throws IOException
     * 
     */
    public void loadGenomeGFFData(String genomeGFFFile) throws IOException{

        
        String line = null;
        BufferedReader brGFF = new BufferedReader(new FileReader(new File(genomeGFFFile)));
            while((line = brGFF.readLine())!= null){
                this.addFeature(line);
            }
        brGFF.close();
        logger.info("read " + genomeFeatureList.size() + " feature entries");
        
    }
    
    
}
