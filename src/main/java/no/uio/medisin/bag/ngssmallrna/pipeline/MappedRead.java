/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

/**
 *
 * @author simonray
 */
public class MappedRead implements Comparable<MappedRead>{
    private int startPos;
    private int endPos;
    private String chr;
    private String strand;
    private int count;

    
    
    public MappedRead(int s, int e, String c, String t, int n){
        startPos = s;
        endPos = e;
        chr = c;
        strand = t;
        count = n;
    }
    
    /**
     * sort by chromosome, strand and then position
     * 
     * @param mappedRead
     * @return 
     */
    @Override
    public int compareTo(MappedRead mappedRead){
        
        if (mappedRead.chr.equals(chr) && mappedRead.strand.equals(strand))
        {
            return startPos - mappedRead.startPos ;
        }
        
        if(mappedRead.strand.equals(strand)==false){
            if(strand.equals("+") && mappedRead.strand.equals("-")){
                return -1;
            }else{
                return +1;
            }            
        }
            
        int l = Math.min(mappedRead.chr.length(), chr.length());
        int i=0;
        while(mappedRead.chr.charAt(i)==(chr.charAt(i)) && i<l){
            i++;
        }
        return  chr.charAt(i) - mappedRead.chr.charAt(i) ;
        
    }
    
    /**
     * print the mapped read information
     * 
     * @return 
     */
    @Override
    public String toString(){
        return    "chr "    + chr + "\t" 
                + "start "  + startPos + "\t"
                + "end "    + endPos + "\t"
                + "strand " + strand + "\t"
                + "count "  + getCount();
    }
    
    
    /**
     * @return the startPos
     */
    public int getStartPos() {
        return startPos;
    }

    /**
     * @return the endPos
     */
    public int getEndPos() {
        return endPos;
    }

    /**
     * @return the chr
     */
    public String getChr() {
        return chr;
    }

    /**
     * @return the strand
     */
    public String getStrand() {
        return strand;
    }

    /**
     * @return the count
     */
    public int getCount() {
        return count;
    }
    
}
