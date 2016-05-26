import os
import re
from l3sdk import define, Process, require

@require(mem_mb=8000,high_io=True  )
class WrapperClass(define.Wrapper):
    class Inputs(define.Inputs):
        #annotation = define.input(name = "Annotation", required = True, description = "Reference sequence")
        #reads = define.input(name = "Aligned  Reads", required = True, description = "Reference sequence")
        reads = define.input(name = 'Input SAM or BAM file',required = True, description = "Input SAM or BAM file."   )

    class Outputs(define.Outputs):
        out = define.output(name = 'File to write the output to.' , description = "Files(s) containing summary alignment metrics."  )

    class Params(define.Params):
        Minimunmappingquality= define.integer( name = "Minimun mapping quality", default = 20, description = "[MINIMUM_MAPPING_QUALITY]Minimum mapping quality for a read to contribute coverage. Default value: 20. This option can be set to 'null' to clear the default value")
        Minimumbasequality= define.integer( name = "Minimum base quality", default = 20, description = "[MINIMUM_BASE_QUALITY]Minimum base quality for a base to contribute coverage. Default value: 20. This option can be set to 'null' to clear the default value.")
        Coveragecap= define.integer( name = "Coverage cap", default = 250, description = "[COVERAGE_CAP]Treat bases with coverage exceeding this value as if they had coverage at this value. Default value: 250. This option can be set to 'null' to clear the default value.")
        Stopafter= define.real( name = "Stop after", default = -1, description = "[STOP_AFTER]For debugging purposes, stop after processing this many genomic bases. Default value: -1. This option can be set to 'null' to clear the default value.")
        Validationstringency= define.enum( name = "Validation stringency", default = "SILENT", values = [( 'STRICT','STRICT','' ),( 'LENIENT','LENIENT','' ),( 'SILENT','SILENT','' )], description = "[VALIDATION_STRINGENCY]Validation stringency for all BAM/SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.")

    def execute(self):
        options=[]        
        if(self.params.Minimunmappingquality):
        #integer#[MINIMUM_MAPPING_QUALITY]Minimum mapping quality for a read to contribute coverage. Default value: 20. This option can be set to 'null' to clear the default value
            options.append( 'MINIMUM_MAPPING_QUALITY='+ str(self.params.Minimunmappingquality))
        if(self.params.Minimumbasequality):
        #integer#[MINIMUM_BASE_QUALITY]Minimum base quality for a base to contribute coverage. Default value: 20. This option can be set to 'null' to clear the default value.
            options.append( 'MINIMUM_BASE_QUALITY='+ str(self.params.Minimumbasequality))
        if(self.params.Coveragecap):
        #integer#[COVERAGE_CAP]Treat bases with coverage exceeding this value as if they had coverage at this value. Default value: 250. This option can be set to 'null' to clear the default value.
            options.append( 'COVERAGE_CAP='+ str(self.params.Coveragecap))
        if(self.params.Stopafter):
        #float#[STOP_AFTER]For debugging purposes, stop after processing this many genomic bases. Default value: -1. This option can be set to 'null' to clear the default value.
            options.append( 'STOP_AFTER='+ str(self.params.Stopafter))
        if(self.params.Validationstringency):
        #enum#Validation stringency for all BAM/SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
            options.append('VALIDATION_STRINGENCY='+self.params.Validationstringency)
          
        #java -jar /opt/bin/picard.jar CollectAlignmentSummaryMetrics INPUT=test-data/ERR315327_.accepted_hits.bam OUTPUT=a.out
        #cmd_run = ['java','-Xmx6g','-jar','/opt/bin/picard.jar','CollectWgsMetrics','REFERENCE_SEQUENCE=/opt/db/human_g1k_v37_decoy.fasta']    
        cmd_run = ['java','-Xmx7g','-jar','/opt/bin/picard.jar','CollectWgsMetrics','REFERENCE_SEQUENCE=/opt/db/human_g1k_v37_decoy.fasta']   

        fileNamePath, fileExtension = os.path.splitext(self.inputs.reads)
        cmd_run2 = []
        cmd_run2.extend(cmd_run)
        cmd_run2.extend(['INPUT='+self.inputs.reads,'OUTPUT='+fileNamePath+'.wgs_metrics.txt'])
        Process(*cmd_run2).run()
            
        self.outputs.out = fileNamePath+'.wgs_metrics.txt'
        self.outputs.out.meta = self.inputs.reads.make_metadata()
           
def test_wrapperclass():
        inputs = {'reads': '/l3bioinfo/test-data/output.mergedsorted.bam' }
        params = {}
        outputs = WrapperClass(inputs, params).test()
        #assert outputs.out.endswith('mock.bam')
