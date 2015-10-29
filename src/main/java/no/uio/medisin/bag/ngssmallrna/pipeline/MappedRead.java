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
public class MappedRead {
    private int startPos;
    private int endPos;
    private String chr;

    
    
    public MappedRead(int s, int e, String c){
        startPos = s;
        endPos = s;
        chr = c;
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
    
}
