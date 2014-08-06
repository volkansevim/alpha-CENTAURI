import os
import sys
pbcore.io.FastaIO import FastqReader, FastaReader
from pbcore.io.BlasrIO import M4Reader
import optparse, logging

class M4FastaFilter:
    def __init__( self ):
        self.__parseArgs( )
        self.__initLog( )

    def __parseArgs( self ):
        """Handle command line argument parsing"""
        
        usage = "%prog [--help] [options] ARG1"
        parser = optparse.OptionParser( usage=usage, description=__doc__ )

        parser.add_option( "-l", "--logFile", help="Specify a file to log to. Defaults to stderr." )
        parser.add_option( "-d", "--debug", action="store_true", help="Increases verbosity of logging" )
        parser.add_option( "-i", "--info", action="store_true", help="Display informative log entries" )
        parser.add_option( "-p", "--profile", action="store_true", help="Profile this script, dumping to <scriptname>.profile" )
        parser.add_option( "-s", "--scoreCut", help="Score cutoff to accept alignment")
        parser.add_option( "--pctCut", help="Percent Similarity cutoff to accept alignment")
        parser.add_option( "--lenCut", help="aligment length cutoff")
        parser.add_option( "--fastq", action="store_true", help="input is fastq" )


        parser.set_defaults( logFile=None, debug=False, info=False, profile=False, scoreCut = 100, pctCut = .7, lenCut = 100, fastq=False)
        
        self.opts, args = parser.parse_args( )

        if len(args) != 1:
            parser.error( "Expected a single argument." )

        self.blasrfn = args[0]
        self.fastafn = args[1]
        self.minScore = int(self.opts.scoreCut) 
        self.minPctId = float(self.opts.minPctCut)
        self.minAlignLen = int(self.opts.lenCut)



    def __initLog( self ):
        """Sets up logging based on command line arguments. Allows for three levels of logging:
        logging.error( ): always emitted
        logging.info( ): emitted with --info or --debug
        logging.debug( ): only with --debug"""

        logLevel = logging.DEBUG if self.opts.debug else logging.INFO if self.opts.info else logging.ERROR
        logFormat = "%(asctime)s [%(levelname)s] %(message)s"
        if self.opts.logFile != None:
            logging.basicConfig( filename=self.opts.logFile, level=logLevel, format=logFormat )
        else:
            logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
                                                                 
    def run( self ):
        """Executes the body of the script."""
    
        logging.info("Log level set to INFO")
        logging.debug("Log Level set to DEBUG")

        readids = readBlasrIDsToDict ()
        filterFastaByIds (readids)
        return 0

    def readBlasrIDsToDict (self):
    	readids = set()
    	for m4record in M4Reader(self.blasrfn):
    		if m4recore.score > self.minScore and m4record.percentSimilarity > self.minPctId and self.tEnd-self.tStart > self.minAlignLen:
    			readids.add(m4record.qName)
    	return readids

    def filterFastaByIds (self, readids):
    	for fastafn in self.fastafns:
    		if self.opts.fastq:
    			for entry in FastaReader (fastafn):
    				if entry.name in readids:
    					print entry
    		else:
    			for entry in FastaReader (fastafn):
    				if entry.name in readids:
    					print entry



if __name__ == "__main__":
    app = M4FastaFilter()
    if app.opts.profile:
        import cProfile
        cProfile.run( 'app.run()', '%s.profile' % sys.argv[0] )
    sys.exit( app.run() )

