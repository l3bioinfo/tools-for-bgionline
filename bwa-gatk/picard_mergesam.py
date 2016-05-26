import os
import re
from l3sdk import define, Process, require

@require(mem_mb=58000  )
class WrapperClass(define.Wrapper):
    class Inputs(define.Inputs):
        #annotation = define.input(name = "Annotation", required = True, description = "Reference sequence")
        #reads = define.input(name = "Aligned  Reads", required = True, description = "Reference sequence")
        reads = define.input(name = 'Input BAM/SAM files',list = True ,required = True, description = "Input SAM or BAM file."   )

    class Outputs(define.Outputs):
        out = define.output(name = 'Merged BAM/SAM files', description = "Merged BAM/SAM files"   )
        ind = define.output(name = 'Index file', alt_path='/extra' )
    class Params(define.Params):
        Assumesorted= define.boolean( name = "Assume sorted", default = True, description = "[ASSUME_SORTED]If true, assume that the input files are in the same sort order as the requested output sort order, even if their headers say otherwise. [Default: false]")
        Sortorder= define.enum( name = "Sort order", default = "coordinate", values = [( 'unsorted','unsorted','' ),( 'queryname','queryname','' ),( 'coordinate','coordinate','' )], description = "[SORT_ORDER]Desired sort order. [default: coordinate]")
        Validationstringency= define.enum( name = "Validation stringency", default = "SILENT", values = [( 'STRICT','STRICT','' ),( 'LENIENT','LENIENT','' ),( 'SILENT','SILENT','' )], description = "Validation stringency for all BAM/SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.")
        MergeSequenceDictionary= define.boolean( name = "Merge Sequence Dictionary", default = False, description = "[MERGE_SEQUENCE_DICTIONARIES] Merge the sequence dictionaries")
        Compressionlevel= define.integer( name = "Compression level", default = 5, description = "Compression level for all compressed files created (e.g. BAM and GELI)")
        CreateIndex= define.boolean( name = "Create Index", default = True, description = "Whether to create a BAM index when writing a coordinate-sorted BAM file")

    def execute(self):
        options=[]        
        if(self.params.Assumesorted):
        #boolean#[ASSUME_SORTED]If true, assume that the input files are in the same sort order as the requested output sort order, even if their headers say otherwise. [Default: false]
            options.append( 'ASSUME_SORTED=true')
        if(self.params.Sortorder):
        #enum#[SORT_ORDER]Desired sort order. [default: coordinate]
            options.append( 'SORT_ORDER='+ self.params.Sortorder)
        if(self.params.MergeSequenceDictionary):
        #boolean#[MERGE_SEQUENCE_DICTIONARIES] Merge the sequence dictionaries
            options.append( 'MERGE_SEQUENCE_DICTIONARIES=true')
        if(self.params.CreateIndex):
            options.append( 'CREATE_INDEX=true')
        #boolean#Whether to create a BAM index when writing a coordinate-sorted BAM file
      
        if(self.params.Validationstringency):
        #enum#Validation stringency for all BAM/SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
            options.append('VALIDATION_STRINGENCY='+self.params.Validationstringency)

        if(self.params.Compressionlevel):
        #integer#Compression level for all compressed files created (e.g. BAM and GELI)
            options.append('COMPRESSION_LEVEL='+str(self.params.Compressionlevel))
            
        #java -jar /opt/bin/picard.jar CollectAlignmentSummaryMetrics INPUT=test-data/ERR315327_.accepted_hits.bam OUTPUT=a.out
        # 2g
        cmd_run = ['java', '-Xmx56g','-jar','/opt/bin/picard.jar','MergeSamFiles','USE_THREADING=true']    
        cmd_run.extend(options)
        for i in range(len(self.inputs.reads)):
            cmd_run.append('INPUT='+self.inputs.reads[i])
            
        fileNamePath, fileExtension = os.path.splitext(self.inputs.reads[0])

        cmd_run.append( 'OUTPUT='+fileNamePath+'sorted'+ fileExtension )
        Process(*cmd_run).run()
            
        self.outputs.out = fileNamePath+'sorted'+ fileExtension  
        self.outputs.out.meta = self.inputs.reads[0].make_metadata()
        
        if (self.params.CreateIndex and os.path.exists(fileNamePath+'sorted'+ fileExtension +'.bai') ):
            self.outputs.ind = fileNamePath+'sorted'+ fileExtension +'.bai' 
            self.outputs.ind.meta = self.inputs.reads[0].make_metadata()        
            
 
def test_wrapperclass():
        inputs = {'reads': ['/l3bioinfo/test-data/ERR315327_.accepted_hits.bam', '/l3bioinfo/test-data/ERR315384_.accepted_hits.bam'] }
          
        params = {}
        outputs = WrapperClass(inputs, params).test()
        #assert outputs.out.endswith('mock.bam')
