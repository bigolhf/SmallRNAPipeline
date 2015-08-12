/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import no.uio.medisin.bag.jmirpara.SimpleSeq;
/**
 * stores an miRBase entry
 * This needs to be made comparable so the entries stored in a list can be sorted
 * to speed things up
 * 
 * @author sr
 */
public class GenomeFeature {

    private static Logger                       logger                      = LogManager.getLogger();


    
    private String chromosome;
    private String source;
    private String featureName;
    private int    startPos;
    private int    endPos;
    private String score;
    private String strand;
    private String sequence;
    private String attributeString;
    

    
    private static final int seqnameCol     = 0;
    private static final int sourceCol      = 1;
    private static final int featureCol     = 2;
    private static final int startyCol       = 3;
    private static final int endCol         = 4;
    private static final int scoreCol       = 5;
    private static final int strandCol      = 6;
    private static final int attribCol      = 7;
    
    

    
    
    
    /**
     * create entry from GFF string
     * 
     *  Col 1 : seqname - name of the chromosome or scaffold
     *  Col 2 : source - name of the program that generated this feature, or the data source (database or project name)
     *  Col 3 : feature - feature type name, e.g. Gene, Variation, Similarity
     *  Col 4 : start - Start position of the feature, with sequence numbering starting at 1
     *  Col 5 : end - End position of the feature, with sequence numbering starting at 1
     *  Col 6 : score - A floating point value.
     *  Col 7 : strand - defined as + (forward) or - (reverse).
     *  Col 8 : frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon,
     *  Col 9 : attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
     * 
     * @param gffString 
     */
    public GenomeFeature(String gffString){
        
        String tokens[] = gffString.split("\t");
        chromosome = tokens[seqnameCol].trim();
        source = tokens[sourceCol].trim();
        featureName = tokens[featureCol].trim();
        startPos = Integer.parseInt(tokens[startyCol].trim());
        endPos = Integer.parseInt(tokens[endCol].trim());
        score = tokens[scoreCol].trim();
        strand = tokens[strandCol].trim();
        attributeString = tokens[attribCol].trim();
        
    }
    
    
    
    /**
     * checks whether Chromosome strings are the same, while attempting
     * to allow for the presence or absence of a variation on the 'Chr' 
     * prefix
     * 
     * @param queryChr
     * @return 
     */
    public Boolean chromosomeMatch(String queryChr){
        
        return GenomeFeature.removeChromosomePrefix(getChromosome()).equals(GenomeFeature.removeChromosomePrefix(queryChr));
        
    }
    
    
    
    
    /**
     * attempt to remove any prefix of the form 'Chr' from the chromosome string
     * 
     * @param chrString
     * @return 
     */
    public static String removeChromosomePrefix(String chrString){
        
        if(chrString.contains("chr")){
            chrString = chrString.replace("chr", "");
        }
        else{
            if(chrString.contains("CHR")){
                chrString = chrString.replace("CHR", "");
            }
            else{
                if(chrString.contains("Chr")){
                    chrString = chrString.replace("Chr", "");
                }
            }
            
        }
        return chrString;
        
    }
        
    
    
    

    
    
    
    
    /**
     * base equality on the mimatID, which should be unique
     * 
     * @param qObject
     * @return 
     */
    @Override
    public boolean equals(Object qObject){
        if (qObject != null && qObject instanceof GenomeFeature)
        {
            return (this.startPos == ((GenomeFeature) qObject).startPos);
        }
        return false;

    }
    
    
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + startPos;
        return result;
    }    

    /**
     * @return the chromosome
     */
    public String getChromosome() {
        return chromosome;
    }

    /**
     * @param chromosome the chromosome to set
     */
    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    /**
     * @return the source
     */
    public String getSource() {
        return source;
    }

    /**
     * @param source the source to set
     */
    public void setSource(String source) {
        this.source = source;
    }

    /**
     * @return the featureName
     */
    public String getFeatureName() {
        return featureName;
    }

    /**
     * @param featureName the featureName to set
     */
    public void setFeatureName(String featureName) {
        this.featureName = featureName;
    }

    /**
     * @return the startPos
     */
    public int getStartPos() {
        return startPos;
    }

    /**
     * @param startPos the startPos to set
     */
    public void setStartPos(int startPos) {
        this.startPos = startPos;
    }

    /**
     * @return the endPos
     */
    public int getEndPos() {
        return endPos;
    }

    /**
     * @param endPos the endPos to set
     */
    public void setEndPos(int endPos) {
        this.endPos = endPos;
    }

    /**
     * @return the score
     */
    public String getScore() {
        return score;
    }

    /**
     * @param score the score to set
     */
    public void setScore(String score) {
        this.score = score;
    }

    /**
     * @return the strand
     */
    public String getStrand() {
        return strand;
    }

    /**
     * @param strand the strand to set
     */
    public void setStrand(String strand) {
        this.strand = strand;
    }

    /**
     * @return the sequence
     */
    public String getSequence() {
        return sequence;
    }

    /**
     * @param sequence the sequence to set
     */
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    /**
     * @return the attributeString
     */
    public String getAttributeString() {
        return attributeString;
    }

    /**
     * @param attributeString the attributeString to set
     */
    public void setAttributeString(String attributeString) {
        this.attributeString = attributeString;
    }
    
       
}
