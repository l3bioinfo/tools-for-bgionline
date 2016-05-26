import os
import os.path
from l3sdk import define, Process, require

@require(mem_mb=10000, cpu=require.CPU_ALL,high_io=True)
class WrapperClass(define.Wrapper):
    class Inputs(define.Inputs):
        reads = define.input(name = "FASTQ Reads", list = True,description="A single FASTQ file for single end or two files for paired end experiments", required=True)
        adapter = define.input(name = "Adapters" , description="Alist of adapters to clip in FASTA format")
        
    class Outputs(define.Outputs):
        out_fq = define.output(name="Filtered FASTQ files",list =True, description="FASTQ files processed according to the set parameters")
        out_skip = define.output(name="Skipped Reads",list =True, description="FASTQ files containing reads that were filtered out from the input")
        out_summary = define.output(name="Summary", description="A text file containing summary of the filtration")
    class Params(define.Params):
        Minimumlengthmatch= define.real( name = "Minimum length match", default = 2.2, description = "[-s]: Log scale for adapter minimum-length-match (2.2)")
        Adapteroccurrencethreshold= define.real( name = "Adapter occurrence threshold", default = 0.25, description = "[-t]: % occurrence threshold before adapter clipping (0.25)")
        Mincliplength= define.integer( name = "Min clip length", default = 1, description = "[-m]: Minimum clip length, overrides scaled auto (1)")
        Maxadapterdifference= define.integer( name = "Max adapter difference", default = 10, description = "[-p]: Maximum adapter difference percentage (10)")
        Setalldefaultparameterstozerodonothing= define.boolean( name = "Set all default parameters to zero/do nothing", default = False, description = "[-0]: default False")
        Minremainingsequencelength= define.integer( name = "Min remaining sequence length", default = 19, description = "[-l]: Minimum remaining sequence length (19)")
        Maxremainingsequencelength= define.integer( name = "Max remaining sequence length", description = "[-L]: Maximum remaining sequence length")
        Removeduplicatereads= define.integer( name = "Remove duplicate reads", default = 0, description = "[-D]: Read_1 has an identical N bases (0)")
        sKewPercentage= define.integer( name = "sKew Percentage", default = 2, description = "[-k]: If any nucleotide is less than the skew percentage, then the whole cycle is removed (2). Set the skew (-k) or N-pct (-x) to 0 to turn it off, this should be done for miRNA, amplicon and other low-complexity situations.")
        Badreadpercentagethreshold= define.integer( name = "Bad read percentage threshold", default = 20, description = "[-x]: 'N' (Bad read) percentage causing cycle removal from ALL read (20). Set the skew (-k) or N-pct (-x) to 0 to turn it off, this should be done for miRNA, amplicon and other low-complexity situations.")
        Qualitythreshold= define.integer( name = "Quality threshold", default = 7, description = "[-q]: Quality threshold causing base removal (7)")
        Trimmingwindowsize= define.integer( name = "Trimming window size", default = 1, description = "[-w]: Window-size for quality trimming (1)")
        Removehomopolymerreads= define.boolean( name = "Remove homopolymer reads", default = False, description = "[-H]: Remove >95% homopolymer reads")
        IlluminaPF= define.boolean( name = "Illumina PF", default = False, description = "[-U|u]: Force disable/enable Illumina PF filtering. Values are -u: disable (default); -U: enable")
        DonttrimNs= define.boolean( name = "Don't trim N's", default = False, description = "[-R]: Don't remove N's from the fronts/ends of reads")
        Subsampling= define.integer( name = "Subsampling", default = 300000, description = "[-C]: Number of reads to use for subsampling (300k)")
        Phredscale= define.integer( name = "Phred-scale", description = "[-P]: Phred-scale (auto-determined)")
        Dontclip= define.boolean( name = "Don't clip", default = False, description = "[-n]: Just output what would be done")
        Onlykeepclippedreads= define.boolean( name = "Only keep clipped reads", default = False, description = "[-K]: Only keep clipped reads")
        Saveskippedreads= define.boolean( name = "Save skipped reads", default = False, description = "[-S]: Output FASTQ files skipped reads on the 'Skipped Reads' output.")
        Minimummeanqualityscore= define.real( name = "Minimum mean quality score", description = "[--qual-mean]: Evaluated after clipping/trimming")
        Minimummeanqualityscoreappliestosecondnonbarcodereadonly= define.real( name = "Minimum mean quality score, applies to second non-barcode read only", description = "[--mate-qual-mean]: Evaluated after clipping/trimming")
        Qualitygreaterthanthreshold= define.string( name = "Quality greater than threshold", description = "[--qual-gt NUM,THR]: Evaluated after clipping/trimming, At least NUM quals > THR")
        Qualitygreaterthanthresholdappliestosecondnonbarcodereadonly= define.string( name = "Quality greater than threshold, applies to second non-barcode read only", description = "[--mate-qual-gt NUM,THR]:Evaluated after clipping/trimming, At least NUM quals > THR")
        MaximumNcallsinareadcanbea= define.real( name = "Maximum N-calls in a read (can be a %)", description = "[--max-ns]: Evaluated after clipping/trimming")
        MaximumNcallsinareadcanbeaappliestosecondnonbarcodereadonly= define.real( name = "Maximum N-calls in a read (can be a %), applies to second non-barcode read only", description = "[--mate-max-ns]: Evaluated after clipping/trimming")
        Homopolymerfilterpercentageasnumber= define.integer( name = "Homopolymer filter percentage, as number", description = "[--homopolymer-pct]: Homopolymer filter percentage, evaluated after clipping/trimming")
        Complexityfilterpercent= define.integer( name = "Complexity filter percent", description = "[--lowcomplex-pct]: Complexity filter percent (95)")
        AdjustcycleCYCnegativeoffsetfromendbyamountAMT= define.string( name = "Adjust cycle CYC (negative - offset from end) by amount AMT", description = "[--cycle-adjust CYC,AMT] Adjust cycle CYC (negative - offset from end) by amount AMT")
        AdjustscoreSCOREbyamountAMT= define.string( name = "Adjust score SCORE by amount AMT", description = "[--phred-adjust SCORE,AMT]: Adjust score SCORE by amount AMT")
        #AutoAdjustToSanger = define.boolean( name = "Auto Adjust to Sanger Scale", default = True, description = "Auto change scale to Sanger")

    def execute(self):
        options=[]
        if(self.params.Minimumlengthmatch):
        #float#[-s]: Log scale for adapter minimum-length-match (2.2)
            options.extend([ '-s', str(self.params.Minimumlengthmatch)])
        if(self.params.Adapteroccurrencethreshold):
        #float#[-t]: % occurrence threshold before adapter clipping (0.25)
            options.extend([ '-t', str(self.params.Adapteroccurrencethreshold)])
        if(self.params.Mincliplength):
        #integer#[-m]: Minimum clip length, overrides scaled auto (1)
            options.extend([ '-m', self.params.Mincliplength])
        if(self.params.Maxadapterdifference):
        #integer#[-p]: Maximum adapter difference percentage (10)
            options.extend([ '-p', self.params.Maxadapterdifference])
        if(self.params.Setalldefaultparameterstozerodonothing):
        #boolean#[-0]: default False
            options.append( '-0')
        if(self.params.Minremainingsequencelength):
        #integer#Minimum remaining sequence length (19)
            options.extend([ '-l', self.params.Minremainingsequencelength])
        if(self.params.Maxremainingsequencelength):
        #integer#[-L]: Maximum remaining sequence length
            options.extend([ '-L', self.params.Maxremainingsequencelength])
        if(self.params.Removeduplicatereads):
        #integer#[-D]: Read_1 has an identical N bases (0)
            options.extend([ '-D', self.params.Removeduplicatereads])
        if(self.params.sKewPercentage):
        #integer#[-k]: If any nucleotide is less than the skew percentage, then the whole cycle is removed (2). Set the skew (-k) or N-pct (-x) to 0 to turn it off, this should be done for miRNA, amplicon and other low-complexity situations.
            options.extend([ '-k', self.params.sKewPercentage])
        if(self.params.Badreadpercentagethreshold):
        #integer#[-x]: 'N' (Bad read) percentage causing cycle removal from ALL read (20). Set the skew (-k) or N-pct (-x) to 0 to turn it off, this should be done for miRNA, amplicon and other low-complexity situations.
            options.extend([ '-x', self.params.Badreadpercentagethreshold])
        if(self.params.Qualitythreshold):
        #integer#[-q]: Quality threshold causing base removal (7)
            options.extend([ '-q', self.params.Qualitythreshold])
        if(self.params.Trimmingwindowsize):
        #integer#[-w]: Window-size for quality trimming (1)
            options.extend([ '-w', self.params.Trimmingwindowsize])
        if(self.params.Removehomopolymerreads):
        #boolean#[-H]: Remove >95% homopolymer reads
            options.append( '-H')
        if(self.params.IlluminaPF):
        #boolean#[-U|u]: Force disable/enable Illumina PF filtering. Values are -u, disable (default), -U, enable
            options.append( '-U')
        if(self.params.DonttrimNs):
        #boolean#[-R]: Don't remove N's from the fronts/ends of reads
            options.append( '-R')
        if(self.params.Subsampling):
        #integer#[-C]: Number of reads to use for subsampling (300k)
            options.extend([ '-C', self.params.Subsampling])
        if(self.params.Phredscale):
        #integer#[-P]: Phred-scale (auto-determined)
            options.extend([ '-P', self.params.Phredscale])
        if(self.params.Dontclip):
        #boolean#[-n]: Just output what would be done
            options.append( '-n')
        if(self.params.Onlykeepclippedreads):
        #boolean#[-K]: Only keep clipped reads
            options.append( '-K')
        if(self.params.Saveskippedreads):
        #boolean#[-S]: Output FASTQ files skipped reads on the 'Skipped Reads' output.
            options.append( '-S')
        if(self.params.Minimummeanqualityscore):
        #float#[--qual-mean]: Evaluated after clipping/trimming
            options.extend([ '--qual-mean', str(self.params.Minimummeanqualityscore)])
        if(self.params.Minimummeanqualityscoreappliestosecondnonbarcodereadonly):
        #float#[--mate-qual-mean]: Evaluated after clipping/trimming
            options.extend([ '--mate-qual-mean', str(self.params.Minimummeanqualityscoreappliestosecondnonbarcodereadonly)])
        if(self.params.Qualitygreaterthanthreshold):
        #string#[--qual-gt NUM,THR]: Evaluated after clipping/trimming, At least NUM quals > THR
            options.extend([ '--qual-gt', self.params.Qualitygreaterthanthreshold])
        if(self.params.Qualitygreaterthanthresholdappliestosecondnonbarcodereadonly):
        #string#[--mate-qual-gt NUM,THR]:Evaluated after clipping/trimming, At least NUM quals > THR
            options.extend([ '--mate-qual-gt', self.params.Qualitygreaterthanthresholdappliestosecondnonbarcodereadonly])
        if(self.params.MaximumNcallsinareadcanbea):
        #float#[--max-ns]: Evaluated after clipping/trimming
            options.extend([ '--max-ns', str(self.params.MaximumNcallsinareadcanbea)])
        if(self.params.MaximumNcallsinareadcanbeaappliestosecondnonbarcodereadonly):
        #float#[--mate-max-ns]: Evaluated after clipping/trimming
            options.extend([ '--mate-max-ns', str(self.params.MaximumNcallsinareadcanbeaappliestosecondnonbarcodereadonly)])
        if(self.params.Homopolymerfilterpercentageasnumber):
        #integer#[--homopolymer-pct]: Homopolymer filter percentage, evaluated after clipping/trimming
            options.extend([ '--homopolymer-pct', self.params.Homopolymerfilterpercentageasnumber])
        if(self.params.Complexityfilterpercent):
        #integer#[--lowcomplex-pct]: Complexity filter percent (95)
            options.extend([ '--lowcomplex-pct', self.params.Complexityfilterpercent])
        if(self.params.AdjustcycleCYCnegativeoffsetfromendbyamountAMT):
        #string#[--cycle-adjust CYC,AMT] Adjust cycle CYC (negative - offset from end) by amount AMT
            options.extend([ '--cycle-adjust', self.params.AdjustcycleCYCnegativeoffsetfromendbyamountAMT])
        if(self.params.AdjustscoreSCOREbyamountAMT):
        #string#[--phred-adjust SCORE,AMT]: Adjust score SCORE by amount AMT
            options.extend([ '--phred-adjust', self.params.AdjustscoreSCOREbyamountAMT])
        
#        if( self.params.AutoAdjustToSanger ):
#            Process('perl' ,'/opt/bin/fastq_detect.pl',self.inputs.reads[0],'1000' ).run()
#        if (  os.path.exists( 'report.txt' ) ):
#            options.extend([ '--phred-adjust',  '-31'])

        run_cmd = ['/opt/bin/ea-utils/fastq-mcf']
        
        is_phed64 = "Phred+33"
        Process('perl' ,'/opt/bin/fastq_detect.pl',self.inputs.reads[0],'1000' ).run()
        if (  os.path.exists( 'report.txt' ) ):
            is_phed64 = "Phred+64"
        
        inputfiles = []
        for x in self.inputs.reads:
            fileNamePath, fileExtension = os.path.splitext(x)
            inputfiles.append(x)
            options.extend(['-o',fileNamePath+'.clip'+fileExtension])
            self.outputs.out_fq.add_file ( fileNamePath+'.clip'+fileExtension)
            self.outputs.out_fq[-1].meta =  x.make_metadata( _quality_scale=is_phed64 )
            # Phred+64 Phred+33
            if (self.params.Saveskippedreads):
                self.outputs.out_skip.add_file ( fileNamePath+'.clip'+fileExtension+'.skip')
                self.outputs.out_skip[-1].meta =  x.make_metadata( _quality_scale=is_phed64 )

        run_cmd.extend(options)

        if (self.inputs.adapter):
            run_cmd.append(self.inputs.adapter)
        else:
            run_cmd.extend(["-f",'/dev/null'])
            
#        options2 = []
#        for x in self.inputs.reads:
#            fileNamePath, fileExtension = os.path.splitext(x)
#            run_cmd.append(x)
#            options2.extend(['-o',fileNamePath+'.clip'+fileExtension])
#            self.outputs.out_fq.add_file ( fileNamePath+'.clip'+fileExtension)
#            self.outputs.out_fq[-1].meta =  x.make_metadata()
#            if (self.params.Saveskippedreads):
#                self.outputs.out_skip.add_file ( fileNamePath+'.clip'+fileExtension+'.skip')
#                self.outputs.out_skip[-1].meta =  x.make_metadata()

        run_cmd.extend(inputfiles)
        #run_cmd.extend(['>',fileNamePath+'.fastq-mcf_summary.txt'  ,'||','true'])
        Process(*run_cmd,stdout= fileNamePath+'.fastq-mcf_summary.txt' ).run()
        self.outputs.out_summary = fileNamePath+'.fastq-mcf_summary.txt' 
        self.outputs.out_summary.meta = self.inputs.reads[0].make_metadata()
        
        
def test_wrapperclass():
        inputs = {'reads': ['/l3bioinfo/test-data/130629_I321_FCC221BACXX_L5_HUM_1.fq.gz'] }
        params = {}
        outputs = WrapperClass(inputs, params).test()        
