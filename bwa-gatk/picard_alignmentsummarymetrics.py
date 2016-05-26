import os
import re
from l3sdk import define, Process, require

@require(mem_mb=2048,high_io=True  )
class WrapperClass(define.Wrapper):
    class Inputs(define.Inputs):
        #annotation = define.input(name = "Annotation", required = True, description = "Reference sequence")
        #reads = define.input(name = "Aligned  Reads", required = True, description = "Reference sequence")
        reads = define.input(name = 'Input SAM or BAM file',list = True ,required = True, description = "Input SAM or BAM file."  )

    class Outputs(define.Outputs):
        out = define.output(name = 'File to write the output to.',list = True , description = "Files(s) containing summary alignment metrics." )

    class Params(define.Params):
        Maximuminsertsize= define.integer( name = "Maximum insert size", default = 100000, description = "[MAX_INSERT_SIZE]Paired end reads above this insert size will be considered chimeric along with inter-chromosomal pairs. [Default: 100000].")
        Adaptersequence= define.string( name = "Adapter sequence", description = "[ADAPTER_SEQUENCE]List of adapter sequences to use when processing the alignment metrics This option may be specified 0 or more times. Separate each by comma or space.")
        Metricaccumulationlevel= define.string( name = "Metric accumulation level", description = "[METRIC_ACCUMULATION_LEVEL]The level(s) at which to accumulate metrics. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} This option may be specified 0 or more times. Separate each by comma or space.")
        Validationstringency= define.enum( name = "Validation stringency", default = "SILENT", values = [( 'STRICT','STRICT','' ),( 'LENIENT','LENIENT','' ),( 'SILENT','SILENT','' )], description = "[VALIDATION_STRINGENCY]Validation stringency for all BAM/SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.")
        Isbisulfitesequenced= define.boolean( name = "Is bisulfite sequenced", default = False, description = "[IS_BISULFITE_SEQUENCED]Whether the SAM or BAM file consists of bisulfite sequenced reads. [Default: false].")
        Assumesorted= define.boolean( name = "Assume sorted", default = True, description = "[ASSUME_SORTED]If true (default), then the sort order in the header file will be ignored. [Default: true].")
        Compressionlevel= define.integer( name = "Compression level", default = 5, description = "[COMPRESSION_LEVEL]Compression level for all compressed files created (e.g. BAM and GELI)")
        CreateIndex= define.boolean( name = "Create Index", default = True, description = "[CREATE_INDEX]Whether to create a BAM index when writing a coordinate-sorted BAM file")

    def execute(self):
        options=[]        
        if(self.params.Maximuminsertsize):
        #integer#Paired end reads above this insert size will be considered chimeric along with inter-chromosomal pairs. [Default: 100000].
            options.append( 'MAX_INSERT_SIZE=' + str(self.params.Maximuminsertsize))
        if(self.params.Adaptersequence):
        #string#This option may be specified 0 or more times.
            for x in re.split(',| ' , self.params.Adaptersequence):
                options.append('ADAPTER_SEQUENCE='+x)
        if(self.params.Metricaccumulationlevel):
        #enum#The level(s) at which to accumulate metrics. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} This option may be specified 0 or more times.
            for x in re.split(',| ' , self.params.Metricaccumulationlevel):
                options.append('METRIC_ACCUMULATION_LEVEL='+x)
        
        if(self.params.Validationstringency):
        #enum#Validation stringency for all BAM/SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
            options.append('VALIDATION_STRINGENCY='+self.params.Validationstringency)

        if(self.params.Isbisulfitesequenced):
        #boolean#Whether the SAM or BAM file consists of bisulfite sequenced reads. [Default: false].
            options.append( 'IS_BISULFITE_SEQUENCED=true')
        if(self.params.Assumesorted==False):
        #boolean#If true (default), then the sort order in the header file will be ignored. [Default: true].
        #default is true
            options.append( 'ASSUME_SORTED=false')

        if(self.params.Compressionlevel):
        #integer#Compression level for all compressed files created (e.g. BAM and GELI)
            options.append('COMPRESSION_LEVEL='+str(self.params.Compressionlevel))
        if(self.params.CreateIndex):
            options.append('CREATE_INDEX=true')
        #boolean#Whether to create a BAM index when writing a coordinate-sorted BAM file
            
        #java -jar /opt/bin/picard.jar CollectAlignmentSummaryMetrics INPUT=test-data/ERR315327_.accepted_hits.bam OUTPUT=a.out
        #cmd_run = ['java','-Xmx750M','-jar','/opt/bin/picard.jar','CollectAlignmentSummaryMetrics']   
        cmd_run = ['java','-Xmx2g','-jar','/opt/bin/picard.jar','CollectAlignmentSummaryMetrics']   
        for i in range(len(self.inputs.reads)):
            fileNamePath, fileExtension = os.path.splitext(self.inputs.reads[i])
            cmd_run2 = []
            cmd_run2.extend(cmd_run)
            cmd_run2.extend(['INPUT='+self.inputs.reads[i],'OUTPUT='+fileNamePath+'.summary_metrics.txt'])
            Process(*cmd_run2).run()
            
            self.outputs.out.add_file( fileNamePath+'.summary_metrics.txt')
            self.outputs.out[-1].meta = self.inputs.reads[i].make_metadata()
           
def test_wrapperclass():
        inputs = {'reads': ['/l3bioinfo/test-data/output.mergedsorted.bam'] }
        params = {}
        outputs = WrapperClass(inputs, params).test()
        #assert outputs.out.endswith('mock.bam')
