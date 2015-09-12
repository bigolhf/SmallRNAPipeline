/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;
/**
 *
 * @author sr
 */
public class RScriptConsole implements RMainLoopCallbacks{
    
    static Logger   logger  = LogManager.getLogger();
    
    @Override
    public void rWriteConsole(Rengine re, String text, int oType) {
        logger.info(text);
    }    

    @Override
    public void   rLoadHistory  (Rengine re, String filename) {
    }			
    
    @Override
    public void   rSaveHistory  (Rengine re, String filename) {
    }			
    
    @Override
    public void   rFlushConsole (Rengine re) {
    }
	
    @Override
    public String rChooseFile(Rengine re, int newFile) {
        return null;        
    }    

    @Override
    public String rReadConsole(Rengine re, String prompt, int addToHistory) {
        return null;
    }

    @Override
    public void rBusy(Rengine re, int which) {
        logger.info("rBusy(" + which + ")");
    }
    
    @Override
    public void rShowMessage(Rengine re, String message) {
        logger.warn("rShowMessage <" + message + ">");
    }
	
    
    public String runScript(){
        
	if (!Rengine.versionCheck()) {
	    logger.info("** Version mismatch - Java files don't match library version.");
	    System.exit(1);
	}
        logger.info("Creating Rengine ...");
        Rengine rEngine=new Rengine (new String [] {"--vanilla"}, false, null);
        logger.info("Rengine created, waiting for R");
		// the engine creates R is a new thread, so we should wait until it's ready
        if (!rEngine.waitForR()) {
            System.out.println("Cannot load R");
            return "fail";
        }
        
      String rScript = "/home/sr/programming/R/helloWorld.R";
      rEngine.eval(String.format("source('%s')", rScript));
      
      
      REXP result = rEngine.eval("greeting    ");
      System.out.println("Greeting R: " + result.asString());
      
      rEngine.end();
      return "success";
    }
    
    public static void main(String[] args){
        
        String s = "";
        RScriptConsole rscriptConsole = new RScriptConsole();
        rscriptConsole.runScript();

    }
}
