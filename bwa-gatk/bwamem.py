import os
import re
from l3sdk import define, Process, require

@require(mem_mb=58000,high_io=True,disk_space_gb=320)
class WrapperClass(define.Wrapper):
    class Inputs(define.Inputs):
        #annotation = define.input(name = "Annotation", required = True, description = "Reference sequence")
        #reads = define.input(name = "Aligned  Reads", required = True, description = "Reference sequence")
        reads = define.input(name = 'Read Sequence',list = True ,required = True, description = "Read Sequence")

    class Outputs(define.Outputs):
        out = define.output(name = 'Aligned SAM or BAM file', description = "SAM or BAM Files(s) after sorting")
        out_bai = define.output(name = 'BAM index file', description = "Bai file")
        
    class Params(define.Params):
        Minimumseedlength= define.integer( name = "Minimum seed length", default = 19, description = "[-k]Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. [default: 19]")
        Bandwidthforbandedalignment= define.integer( name = "Band width for banded alignment", default = 100, description = "[-w]Band width in the banded alignment [default: 100]")
        OffdiagonalXdropoff= define.integer( name = "Off-diagonal X-dropoff", default = 100, description = "[-d]Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. [default: 100]")
        TriggerreseedingforaMEMlongerthanminSeedLenFLOAT= define.real( name = "Trigger re-seeding for a MEM longer than minSeedLen*FLOAT", default = 1.5, description = "[-r]This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. [default: 1.5]")
        SkipseedswithmorethanINToccurrences= define.integer( name = "Skip seeds with more than INT occurrences", default = 500, description = "[-c]Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. [default: 500]")
        Dropchainfraction= define.real( name = "Drop chain fraction", default = 0.5, description = "[-D]Drop chains shorter than FLOAT fraction of the longest overlapping chain.")
        Dropchainlength= define.integer( name = "Drop chain length", default = 0, description = "[-W]Discard a chain if seeded bases shorter than INT.")
        Materescuerounds= define.integer( name = "Mate rescue rounds", default = 50, description = "[-m] Perform at most INT rounds of mate rescues for each read.")
        Skipmaterescue= define.boolean( name = "Skip mate rescue", default = False, description = "[-S] Skip mate rescue")
        SkippairingmaterescueperformedunlessSalsoinuse= define.boolean( name = "Skip pairing; mate rescue performed unless -S also in use", default = False, description = "[-P] In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.")
        Discardexactmatches= define.boolean( name = "Discard exact matches", default = False, description = "[-e] Discard full-length exact matches")
        Readtype= define.enum( name = "Read type", default = "None", values = [( 'None','None','' ),( 'pacbio','pacbio','' ),( 'pbread','pbread','' )], description = "[-x] Read type. Setting -x changes multiple parameters unless overridden pacbio: -k17 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -L0; pbread: -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -N25 -FeaD.001")
        Scoreforasequencematch= define.integer( name = "Score for a sequence match", default = 1, description = "[-A] Score for a sequence match. [default: 1]")
        Penaltyforamismatch= define.integer( name = "Penalty for a mismatch", default = 4, description = "[-B] Penalty for a mismatch. [default: 4]")
        Gapopenpenaltyfordeletions= define.integer( name = "Gap open penalty for deletions", default = 6, description = "[-O] Gap open penalty for deletions [default: 6]")
        Gapopenpenaltyforinsertions= define.integer( name = "Gap open penalty for insertions", default = 6, description = "[-O] Gap open penalty for insertions [default: 6]")
        Gapextensionpenaltyfordeletion= define.integer( name = "Gap extension penalty for deletion", default = 1, description = "[-E] Gap extension penalty for deletion. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). [default: 1]")
        Gapextensionpenaltyforinsertion= define.integer( name = "Gap extension penalty for insertion", default = 1, description = "[-E] Gap extension penalty for insertion. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). [default: 1]")
        Penaltyfor5endclipping= define.integer( name = "Penalty for 5'-end clipping", default = 5, description = "[-L] When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. [default: 5]")
        Penaltyfor3endclipping= define.integer( name = "Penalty for 3'-end clipping", default = 5, description = "[-L] When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. [default: 5]")
        Penaltyforanunpairedreadpair= define.integer( name = "Penalty for an unpaired read pair", default = 17, description = "[-U] BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing. [default: 17]")
        Firstqueryfileconsistsofinterleavedpairedendsequences= define.boolean( name = "First query file consists of interleaved paired-end sequences", default = False, description = "Assume the first input query file is interleaved paired-end FASTA/Q. See the command description for details.")
        XAtag= define.integer( name = "XA tag", default = 5, description = "[-h]If #hits < INT, output all in the XA tag")
        Scorethreshold= define.integer( name = "Score threshold", default = 30, description = "[-T]Minimum score to output [default: 30]")
        OutputallalignmentsforSEorunpairedPE= define.boolean( name = "Output all alignments for SE or unpaired PE", default = False, description = "[-a]Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.")
        AppendappendFASTAQcommenttoSAMoutput= define.boolean( name = "Append append FASTA/Q comment to SAM output", default = False, description = "[-C]This option can be used to transfer read meta information (e.g. barcode) to the SAM output. Note that the FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead to incorrect SAM output.")
        Usesoftclippingforsupplementaryalignments= define.boolean( name = "Use soft clipping for supplementary alignments.", default = False, description = "[-Y]Use soft clipping for supplementary alignments.")
        Markshortersplithitsassecondary= define.boolean( name = "Mark shorter split hits as secondary.", default = True, description = "[-M]Mark shorter split hits as secondary (for Picard compatibility).")
        Completereadgroupheaderline= define.string( name = "Complete read group header line", description = "[-R]Specify the read group in a format like '@RG\tID:foo\tSM:bar'. This value takes precedence over per-attribute parameters. [default: constructed from per-attribute parameters or inferred from metadata]")
        
        Outputformat= define.enum( name = "Output format", default = "Sorted BAM", values = [( 'SAM','SAM','' ),( 'BAM','BAM','' ),( 'Sorted BAM','Sorted BAM','' )], description = "Select format to output. Sorted BAM option will output coordinate sorted BAM.")
        CreateIndex = define.boolean(name = "Create Index",default = True,description="Create Index for Sorted BAM file")
        Filteroutsecondaryalignments= define.boolean( name = "Filter out secondary alignments", default = True, description = "Set to true to filter out secondary alignments. Works only with output format set to BAM or Sorted BAM")
        Duplication= define.enum( name = "Duplication", default = "None", values = [( 'None','None','' ),( 'Mark Duplicates','Mark Duplicates','' ),( 'Remove duplicates','Remove duplicates','' )], description = "Only works for Sorted BAM output. Remove duplicates reads from all output files. Implies: Exclude reads marked as duplicates from discordant, splitter, and/or unmapped file.")
        #SorternumberofGBs= define.integer( name = "Sorter - number of GBs", default = 0, description = "If set to zero, auto-detect best algorithm, else set desired value. [default: 0]")
        #SplitfileslargerthanGB= define.integer( name = "Split files larger than GB", default = 150, description = "Files larger than this value will be split, into this sized chunks for alignment.This number is considered for compressed (.gz) files. For uncompressed files a 3x larger value will be taken. This value is in GB.")

    def execute(self):
        options=[]
        if(self.params.Minimumseedlength):
        #integer#[-k]Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. [default: 19]
            options.extend([ '-k', self.params.Minimumseedlength])
        if(self.params.Bandwidthforbandedalignment):
        #integer#[-w]Band width in the banded alignment [default: 100]
            options.extend([ '-w', self.params.Bandwidthforbandedalignment])
        if(self.params.OffdiagonalXdropoff):
        #integer#[-d]Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. [default: 100]
            options.extend([ '-d', self.params.OffdiagonalXdropoff])
        if(self.params.TriggerreseedingforaMEMlongerthanminSeedLenFLOAT):
        #float#[-r]This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. [default: 1.5]
            options.extend([ '-r', str(self.params.TriggerreseedingforaMEMlongerthanminSeedLenFLOAT)])
        if(self.params.SkipseedswithmorethanINToccurrences):
        #integer#[-c]Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. [default: 500]
            options.extend([ '-c', self.params.SkipseedswithmorethanINToccurrences])
        if(self.params.Dropchainfraction):
        #float#[-D]Drop chains shorter than FLOAT fraction of the longest overlapping chain.
            options.extend([ '-D', str(self.params.Dropchainfraction)])
        if(self.params.Dropchainlength):
        #integer#[-W]Discard a chain if seeded bases shorter than INT.
            options.extend([ '-W', self.params.Dropchainlength])
        if(self.params.Materescuerounds):
        #integer#[-m] Perform at most INT rounds of mate rescues for each read.
            options.extend([ '-m', self.params.Materescuerounds])
        if(self.params.Skipmaterescue):
        #boolean#[-S] Skip mate rescue
            options.append( '-S')
        if(self.params.SkippairingmaterescueperformedunlessSalsoinuse):
        #boolean#[-P] In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.
            options.append( '-P')
        if(self.params.Discardexactmatches):
        #boolean#[-e] Discard full-length exact matches
            options.append( '-e')
        if(self.params.Readtype and self.params.Readtype!="None"):
        #enum#[-x] Read type. Setting -x changes multiple parameters unless overridden pacbio: -k17 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -L0; pbread: -k13 -W40 -c1000 -r10 -A2 -B5 -O2 -E1 -N25 -FeaD.001
            options.extend([ '-x', self.params.Readtype])
        if(self.params.Scoreforasequencematch):
        #integer#[-A] Score for a sequence match. [default: 1]
            options.extend([ '-A', self.params.Scoreforasequencematch])
        if(self.params.Penaltyforamismatch):
        #integer#[-B] Penalty for a mismatch. [default: 4]
            options.extend([ '-B', self.params.Penaltyforamismatch])
        
        if(self.params.Gapopenpenaltyfordeletions and self.params.Gapopenpenaltyforinsertions):
        #integer#[-O] Gap open penalty for deletions [default: 6]
            options.extend([ '-O', str(self.params.Gapopenpenaltyfordeletions)+","+ str(self.params.Gapopenpenaltyforinsertions)])
        else:
            if(self.params.Gapopenpenaltyfordeletions ):
                #integer#[-O] Gap open penalty for deletions [default: 6]
                options.extend([ '-O', self.params.Gapopenpenaltyfordeletions])
            if(self.params.Gapopenpenaltyforinsertions):
                #integer#[-O] Gap open penalty for insertions [default: 6]
                options.extend([ '-O', self.params.Gapopenpenaltyforinsertions])
                
        if(self.params.Gapextensionpenaltyfordeletion and self.params.Gapextensionpenaltyforinsertion):
        #integer#[-O] Gap open penalty for deletions [default: 6]
            options.extend([ '-E', str(self.params.Gapextensionpenaltyfordeletion)+","+ str(self.params.Gapextensionpenaltyforinsertion)])
        else:
            if(self.params.Gapextensionpenaltyfordeletion):
                #integer#[-E] Gap extension penalty for deletion. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). [default: 1]
                options.extend([ '-E', self.params.Gapextensionpenaltyfordeletion])
            if(self.params.Gapextensionpenaltyforinsertion):
                #integer#[-E] Gap extension penalty for insertion. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). [default: 1]
                options.extend([ '-E', self.params.Gapextensionpenaltyforinsertion])

        if(self.params.Penaltyfor5endclipping and self.params.Penaltyfor3endclipping):
        #integer#[-O] Gap open penalty for deletions [default: 6]
            options.extend([ '-L', str(self.params.Penaltyfor5endclipping)+","+ str(self.params.Penaltyfor3endclipping)])
        else:       
            if(self.params.Penaltyfor5endclipping):
                #integer#[-L] When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. [default: 5]
                options.extend([ '-L', self.params.Penaltyfor5endclipping])
            if(self.params.Penaltyfor3endclipping):
                #integer#[-L] When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. [default: 5]
                options.extend([ '-L', self.params.Penaltyfor3endclipping])
        
        if(self.params.Penaltyforanunpairedreadpair):
        #integer#[-U] BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing. [default: 17]
            options.extend([ '-U', self.params.Penaltyforanunpairedreadpair])
#        if(self.params.Firstqueryfileconsistsofinterleavedpairedendsequences):
        #boolean#Assume the first input query file is interleaved paired-end FASTA/Q. See the command description for details.
        if(self.params.XAtag):
        #integer#[-h]If #hits < INT, output all in the XA tag
            options.extend([ '-h', self.params.XAtag])
        if(self.params.Scorethreshold):
        #integer#[-T]Minimum score to output [default: 30]
            options.extend([ '-T', self.params.Scorethreshold])
        if(self.params.OutputallalignmentsforSEorunpairedPE):
        #boolean#[-a]Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
            options.append( '-a')
        if(self.params.AppendappendFASTAQcommenttoSAMoutput):
        #boolean#[-C]This option can be used to transfer read meta information (e.g. barcode) to the SAM output. Note that the FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead to incorrect SAM output.
            options.append( '-C')
        if(self.params.Usesoftclippingforsupplementaryalignments):
        #boolean#[-Y]Use soft clipping for supplementary alignments.
            options.append( '-Y')
        if(self.params.Markshortersplithitsassecondary):
        #boolean#[-M]Mark shorter split hits as secondary (for Picard compatibility).
            options.append( '-M')
        if(self.params.Completereadgroupheaderline):
        #string#[-R]Specify the read group in a format like '@RG\tID:foo\tSM:bar'. This value takes precedence over per-attribute parameters. [default: constructed from per-attribute parameters or inferred from metadata]
            options.extend([ '-R', self.params.Completereadgroupheaderline])
        else:
            if 'ReadGroup' in self.inputs.reads[0].meta:
                sampleName = self.inputs.reads[0].meta.get('ReadGroup')
            elif 'Readgroup' in self.inputs.reads[0].meta:
                sampleName = self.inputs.reads[0].meta.get('Readgroup')
            elif 'readgroup' in self.inputs.reads[0].meta:
                sampleName = self.inputs.reads[0].meta.get('readgroup')
            elif 'RG' in self.inputs.reads[0].meta:
                sampleName = self.inputs.reads[0].meta.get('RG')
            elif 'rg' in self.inputs.reads[0].meta:
                sampleName = self.inputs.reads[0].meta.get('rg')
            elif 'Rg' in self.inputs.reads[0].meta:
                sampleName = self.inputs.reads[0].meta.get('Rg')
            elif 'rG' in self.inputs.reads[0].meta:
                sampleName = self.inputs.reads[0].meta.get('rG')
            else:
                sampleName =  os.path.splitext(os.path.basename(self.inputs.reads[0]))[0]

            if 'SampleName' in self.inputs.reads[0].meta:
                smN = self.inputs.reads[0].meta.get('SampleName')
            elif 'sampleName' in self.inputs.reads[0].meta:
                smN = self.inputs.reads[0].meta.get('sampleName')
            elif 'Samplename' in self.inputs.reads[0].meta:
                smN = self.inputs.reads[0].meta.get('Samplename')
            elif 'samplename' in self.inputs.reads[0].meta:
                smN = self.inputs.reads[0].meta.get('samplename')
            elif 'sample' in self.inputs.reads[0].meta:
                smN = self.inputs.reads[0].meta.get('sample')
            elif 'Sample' in self.inputs.reads[0].meta:
                smN = self.inputs.reads[0].meta.get('Sample')
            else:
                smN = 'DefaultSampleName'
            options.extend(  ['-R', '\"@RG\tID:'+ re.split(' ',sampleName)[0]  + '\tSM:' + re.split(' ',smN)[0] + '\tPL:ILLUMINA\"' ] )


#        if(self.params.Outputformat):
            
        #enum#Select format to output. Sorted BAM option will output coordinate sorted BAM.
#        if(self.params.Filteroutsecondaryalignments):
        #boolean#Set to true to filter out secondary alignments. Works only with output format set to BAM or Sorted BAM
#        if(self.params.Duplication):
        #enum#Remove duplicates reads from all output files. Implies: Exclude reads marked as duplicates from discordant, splitter, and/or unmapped file.
#        if(self.params.SorternumberofGBs):
        #integer#If set to zero, auto-detect best algorithm, else set desired value. [default: 0]
#            options.extend([ 'default: 0', self.params.SorternumberofGBs])
#        if(self.params.SplitfileslargerthanGB):
        #integer#Files larger than this value will be split, into this sized chunks for alignment.This number is considered for compressed (.gz) files. For uncompressed files a 3x larger value will be taken. This value is in GB.
        
        #bwa mem [options] <idxbase> <in1.reads>
        run_cmd = ["/opt/bin/bwa","mem",'-t','12']
        run_cmd.extend(options)
        run_cmd.extend(['/opt/db/human_g1k_v37_decoy.fasta'])
        for i in range(len(self.inputs.reads)):
            run_cmd.append( self.inputs.reads[i] )
            
        fileNamePath, fileExtension = os.path.splitext(self.inputs.reads[0])
        #run_cmd.extend(['>',fileNamePath+'.sam'])
        # 
        if ( self.params.Outputformat != 'BAM' and self.params.Outputformat != 'Sorted BAM'):
            Process(*run_cmd , stdout=(fileNamePath+'.sam') ).run()
        
        filter1 = []
        
        if(self.params.Outputformat=='BAM' or self.params.Outputformat=='Sorted BAM'):
            options2 = ['|', '/opt/bin/bfr', '-b', '256M' ]
            options2.extend ( ['|', '/opt/bin/sambamba_v0.4.7','view','--sam-input','-f','bam','-t','2', ] )
            if(self.params.Filteroutsecondaryalignments):
                filter1.append("not secondary_alignment")
            if (self.params.Duplication):
                filter1.append("not duplicate")
            if (self.params.Filteroutsecondaryalignments or self.params.Duplication):
                options2.extend(['-F', '\"' + ' and '.join(filter1) + '\"'])
            
            options2.extend(['-o',fileNamePath+'.bam','/dev/stdin'])
            run_cmd.extend( options2 )
            Process('echo', '#!/bin/bash', stdout='run_haha.sh').run()
            Process('echo', *run_cmd, stdout='run_bwa.sh').run()
            Process('echo', 'wait', stdout='run_end.sh').run()
            Process('cat', 'run_haha.sh', 'run_bwa.sh', 'run_end.sh', stdout='/l3bioinfo/run.sh').run()
            Process('chmod', '777', '/l3bioinfo/run.sh').run()
            Process('/l3bioinfo/run.sh').run()
            #Process('rm','-f',fileNamePath+'.sam' ).run()

            if (self.params.Outputformat=='BAM'):
                self.outputs.out = fileNamePath+'.bam'
                self.outputs.out.meta = self.inputs.reads[0].make_metadata() 
            else:
                Process('/opt/bin/sambamba_v0.4.7','sort','-m','50G','-t','12', '--tmpdir='+'./temp','-o',fileNamePath+'.sorted.bam',fileNamePath+'.bam').run()
                Process('rm','-f',fileNamePath+'.bam' ).run()

                if(self.params.Duplication=='Mark Duplicates' or self.params.Duplication== 'Remove duplicates' ):
                    temp_options = ['/opt/bin/sambamba_v0.4.7', 'markdup','-t','12']
                    if (self.params.Duplication == 'Remove duplicates'):
                        temp_options.append('--remove-duplicates')
                        temp_options.extend([ fileNamePath+'.sorted.bam', fileNamePath+'.bam' ])
                    Process( *temp_options ).run()
                    self.outputs.out = fileNamePath+'.bam'
                    self.outputs.out.meta = self.inputs.reads[0].make_metadata() 
                else:
                    Process('mv',fileNamePath+'.sorted.bam',fileNamePath+'.bam' ).run()
                    self.outputs.out = fileNamePath+'.bam'
                    self.outputs.out.meta = self.inputs.reads[0].make_metadata() 
                    
                if(self.params.CreateIndex):
                    Process('/opt/bin/sambamba_v0.4.7', 'index','-t','12',fileNamePath+'.bam',fileNamePath+'.bai').run()
                    self.outputs.out_bai = fileNamePath+'.bai'
        else:
            self.outputs.out = fileNamePath+'.sam'
            self.outputs.out.meta = self.inputs.reads[0].make_metadata()
        #sambamba view --sam-input -F "not secondary_alignment" -f bam -o sambamba_out.bam inputfilename
        #sambamba sort  -o sort.bam inputfilename
        

        #if(self.params.Outputformat!='SAM'):
        #    Process('samtools','view').run()
def test_wrapperclass():
        #inputs = {'reads': ['/l3bioinfo/test-data/ERR315384_1.fastq.clip.gz','/l3bioinfo/test-data/ERR315384_2.fastq.clip.gz'] }
        inputs = {'reads': ['/l3bioinfo/test-data/test_150bp_100_1.fq', '/l3bioinfo/test-data/test_150bp_100_2.fq']}
        #inputs = {'reads': ['/l3bioinfo/test-data/20507_I295_FCD109WACXX_L3_SZAXPI008600_1.fq.gz','/l3bioinfo/test-data/20507_I295_FCD109WACXX_L3_SZAXPI008600_2.fq.gz'] }
        #inputs = {'reads': ['/l3bioinfo/test-data/ER315327_1.fastq.gz','/l3bioinfo/test-data/ER315327_2.fastq.gz'] }
        params = {'Duplication':'Remove duplicates', 'Outputformat': 'Sorted BAM'}
        outputBAMs = WrapperClass(inputs, params).test()
        #assert outputs.out.endswith('mock.bam')
