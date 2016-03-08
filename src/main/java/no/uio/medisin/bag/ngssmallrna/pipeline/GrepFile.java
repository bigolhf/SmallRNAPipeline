/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import org.apache.logging.log4j.LogManager;
import static org.apache.logging.log4j.LogManager.getLogger;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author simonray
 */
public class GrepFile {

    private static Logger               logger                      = LogManager.getLogger();
    
    
    public static ArrayList<String> grepFile(String searchString, String fileName) throws IOException{

        ArrayList<String>  cmdStdErr;
        try{
            
            String cmdGrep = "grep + '" + searchString + "' " + fileName;
        
            Runtime rt = Runtime.getRuntime();
            Process proc = rt.exec(cmdGrep);
            BufferedReader brAStdin = new BufferedReader(new InputStreamReader(proc.getInputStream()));
            BufferedReader brAStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));

            String line = null;
            logger.info("<OUTPUT>");
            while ((line = brAStdin.readLine()) != null) {
                logger.info(line);
            }
            logger.info("</OUTPUT>");

            logger.info("<ERROR>");
            int skipCount = 0;


            cmdStdErr = new ArrayList<>();
            while ((line = brAStdErr.readLine()) != null) {
                logger.info(line);
            }

            int exitVal = proc.waitFor();
            logger.info("Process exitValue: " + exitVal);

            brAStdin.close();
            brAStdErr.close();
        }
        catch(IOException | InterruptedException ex ){
            getLogger().info("error executing bsmap command:\n" + ex.toString());
            getLogger().error("error executing bsmap command:\n" + ex.toString());
            throw new IOException("error executing bsmap command");
        }
        
        return cmdStdErr;
    }
}
