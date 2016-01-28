/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

/**
 * Handles processing of a string that contains strand information
 * this is complicated by the different terminology that is used
 * e.g. 5'/3', 3/5, 3p/5p, +/-
 * this class attempts to parse a string and guess the strand
 * 
 * @author simonray
 */
public class StrandString {
    
    public static final     String  PLUSSTRAND      = "+";
    public static final     String  MINUSSTRAND     = "-";
    public static final     String  UNKSTRAND       = "?";
    
    public static Strand guessStrand(String strandString){
        
        if (strandString.contains(PLUSSTRAND)) return Strand.PLUS;
        if (strandString.contains(MINUSSTRAND))  return Strand.MINUS;
        
        return Strand.UNKNOWN;
    }
    
}
