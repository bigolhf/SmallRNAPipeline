/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

/**
 * stores an miRBase entry
 * This needs to be made comparable so the entries stored in a list can be sorted
 * to speed things up
 * 
 * @author sr
 */
public class MiRNAFeature {
    

    private String mimatID;
    private String name;
    private String parent;
    private String chromosome;
    private int    startPos;
    private int    endPos;
    private String strand;
    private String sequence;

    
    /**
     * Constructor specifying name, location and sequence of an miRNA
     * (ideal for predicted miRNAs)
     * 
     * @param n String      name
     * @param c String      chromosome
     * @param s int         start position
     * @param e int         end position
     * @param t String      Strand
     * 
     */
    public MiRNAFeature(String n, String c, int s, int e, String t){
        this(n, c, s, e, t, "", "");
    }
    
    
    
    /**
     * Constructor specifying name, location and sequence of an miRNA
     * (ideal for miRBase entry miRNAs)
     * 
     * @param n String      name
     * @param c String      chromosome
     * @param s int         start position
     * @param e int         end position
     * @param t String      Strand
     * @param m String      MIMAT (miRBase) ID
     * @param p String      MI (miRBase) Parent ID
     * 
     */
    public MiRNAFeature(String n, String c, int s, int e, String t, String m, String p){
        name = n;
        chromosome = c;
        startPos = s;
        endPos = e;
        strand = t;
        parent = p;
        mimatID = m;
        
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
        
        return removeChromosomePrefix(chromosome).equals(removeChromosomePrefix(queryChr));
        
    }
    
    
    
    
    /**
     * attempt to remove any prefix of the form 'Chr' from the chromosome string
     * 
     * @param chrString
     * @return 
     */
    public String removeChromosomePrefix(String chrString){
        
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
     * 
     * @param miRFeat
     * @return
     */
    @Deprecated
    public int compareTo(MiRNAFeature miRFeat) {

        int thisChr = -1;
        int queryChr = -1;
        if(this.getChromosome().contains("chr")){
            thisChr = Integer.parseInt(this.getChromosome().replace("chr", ""));
            queryChr = Integer.parseInt(miRFeat.getChromosome().replace("chr", ""));
        }
        
        if(this.getChromosome().contains("CHR")){
            thisChr = Integer.parseInt(this.getChromosome().replace("CHR", ""));
            queryChr = Integer.parseInt(miRFeat.getChromosome().replace("CHR", ""));
        }
        
        if(this.getChromosome().contains("Chr")){
            thisChr = Integer.parseInt(this.getChromosome().replace("Chr", ""));
            queryChr = Integer.parseInt(miRFeat.getChromosome().replace("Chr", ""));
        }
        
     return thisChr - queryChr;

     }
   
    
    /**
     * @return the mimatID
     */
    public String getMimatID() {
        return mimatID;
    }

    /**
     * @param mimatID the mimatID to set
     */
    public void setMimatID(String mimatID) {
        this.mimatID = mimatID;
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the parent
     */
    public String getParent() {
        return parent;
    }

    /**
     * @param parent the parent to set
     */
    public void setParent(String parent) {
        this.parent = parent;
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
    
}
