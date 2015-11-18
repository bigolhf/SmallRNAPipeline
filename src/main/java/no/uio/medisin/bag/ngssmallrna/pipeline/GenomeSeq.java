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
import java.util.Iterator;
import no.uio.medisin.bag.jmirpara.SimpleSeq;

/**
 *
 * @author simonray
 */
public class GenomeSeq {
    private String genomeName;
    private ArrayList<SimpleSeq> genomeSeq = new ArrayList<>();
    private int noOfBases;
    
    
    public GenomeSeq(String name){
        genomeName = name;
    }
    
    
    
    
    
    
    /**
     * load a genome in FASTA format
     * 
     * @param filename
     * @return
     * @throws IOException 
     */
    public int readFastaGenome(String filename) throws IOException{
        
        
        String line = "";
        String headerLine = "";
        StringBuilder seq = new StringBuilder();
        
        int totalBases = 0;

        BufferedReader brGS = new BufferedReader(new FileReader(new File(filename)));
            headerLine = brGS.readLine();
            if(headerLine.substring(0, 1).equals(">") == false){
                throw new java.io.IOException("bad FASTA format on first line\n");
            }
            headerLine = headerLine.substring(1, headerLine.length());

            while((line = brGS.readLine()) != null){
                totalBases+= line.length()+1;
                if(line.startsWith(">")==true){
                    genomeSeq.add(new SimpleSeq(headerLine, seq.toString()));
                    noOfBases += seq.toString().length();
                    headerLine = line.substring(1, line.length());
                    seq = new StringBuilder();
                }
                else{
                    line=line.replaceAll("\\s+", "").replace('T', 'u').toUpperCase();
                    seq.append(line);
                }
            }
            genomeSeq.add(new SimpleSeq(headerLine, seq.toString()));
            noOfBases += seq.toString().length();
        brGS.close();
        
        return genomeSeq.size();
    }

    
    
    
    /**
     * return subsequence within specified feature
     * 
     * @param featureID
     * @param b
     * @param e
     * @return 
     */
    public String getSubSeq(String featureID, int b, int e){
        Iterator itGN = genomeSeq.iterator();
        while(itGN.hasNext()){
            SimpleSeq simpleSeq = (SimpleSeq) itGN.next();
            if(simpleSeq.getId().equals(featureID))
                return simpleSeq.getSequence(b, e);
        }
        return null;
    }
    
    
    
    
    
    /**
     * @return the genomeName
     */
    public String getGenomeName() {
        return genomeName;
    }
    
    /**
     * @param genomeName the genomeName to set
     */
    public void setGenomeName(String genomeName) {
        this.genomeName = genomeName;
    }
    
    
    
    public static void main(String args[]){
        GenomeSeq genomeSeq = new GenomeSeq("new genome");
        try{
            genomeSeq.readFastaGenome("/data/ngsdata/sweden/test.fa");
            
            System.out.printf(SimpleSeq.complement(genomeSeq.getSubSeq("S2", 5, 17)));
            System.out.printf("done");
        }
        catch(IOException exIO){
            System.err.printf(exIO.toString());
        }
    }

    /**
     * @return the noOfBases
     */
    public int getNoOfBases() {
        return noOfBases;
    }
    
    
    
    
    /**
     * return the number of chromosomes in the genome
     * This may not be chromosomes but rather the number of individual FASTA entries 
     * that were contained in the input file
     * 
     * @return 
     */
    public int getNoOfChr(){
        return genomeSeq.size();
    }
}
