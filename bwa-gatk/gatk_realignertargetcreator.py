import os
import re
from l3sdk import define, Process, require

@require(mem_mb=18000,high_io=True ,disk_space_gb=150)
class WrapperClass(define.Wrapper):
    class Inputs(define.Inputs):
        bam = define.input(name="Read sequence", required=True,description = "Read sequence in BAM format", list=True  ) # , alt_path='/extra'
        bai = define.input(name="BAM index files", required=False,description="According BAM index",list=True )
        exclude_intervals = define.input(name="Exclude Intervals",list=True, required=False   )
        exome_bed = define.input(name="Genomic Intervals",list=True, required=False  )
        Gatk_key =  define.input(name="GATK key", required=False  )

    class Outputs(define.Outputs):
        dedup_bam_intervals = define.output(description="Dedup BAM Intervals", required=True, list=True   ) # , alt_path='/extra'

    class Params(define.Params):
        DisableRandomization= define.boolean( name = "Disable Randomization", default = False, description = "[-ndrs]Completely eliminates randomization from nondeterministic methods. To be used mostly in the testing framework where dynamic parallelism can result in differing numbers of calls to the generator.")
        AllowPotentiallyMisencodedQuals= define.boolean( name = "Allow Potentially Misencoded Quals", default = False, description = "[-allowPotentiallyMisencodedQuals] Do not fail when encountered base qualities that are too high and seemingly indicate a problem with the base quality encoding of the BAM file.")
        BAQCalculationType= define.enum( name = "BAQ Calculation Type", default = "OFF", values = [( 'OFF','OFF','' ),( 'CALCULATE_AS_NECESSARY','CALCULATE_AS_NECESSARY','' ),( 'RECALCULATE','RECALCULATE','' )], description = "[-baq]Type of BAQ calculation to apply in the engine.")
        BAQGapOpenPenalty= define.real( name = "BAQ Gap Open Penalty", default = 40, description = "[-baqGOP]BAQ gap open penalty (Phred Scaled). Default value is 40. 30 is perhaps better for whole genome call sets")
        DefaultBaseQualities= define.integer( name = "Default Base Qualities", default = -1, description = "If reads are missing some or all base quality scores, this value will be used for all base quality scores")
        DisableIndelQuals= define.boolean( name = "Disable Indel Quals", default = False, description = "[-DBQ]If 'true', disables printing of base insertion and base deletion tags (with -BQSR). Turns off printing of the base insertion and base deletion tags when using the -BQSR argument and only the base substitution qualities will be produced.")
        DownsampletoCoverage= define.integer( name = "Downsample to Coverage", description = "[-dcov]Coverage to downsample to at any given locus; note that downsampled reads are randomly selected from all possible reads at a locus. For non-locus-based traversals (eg., ReadWalkers), this sets the maximum number of reads at each alignment start position.")
        DownsampletoFraction= define.real( name = "Downsample to Fraction", description = "[-dfrac]Fraction [0.0-1.0] of reads to downsample to")
        DownsamplingType= define.enum( name = "Downsampling Type", default = "null", values = [( 'NONE','NONE','' ),( 'ALL_READS','ALL_READS','' ),( 'BY_SAMPLE','BY_SAMPLE','' ),('null','null','')], description = "[-dt]Type of reads downsampling to employ at a given locus. Reads will be selected randomly to be removed from the pile based on the method described here")
        EmitOriginalQuals= define.boolean( name = "Emit Original Quals", default = False, description = "[-EOQ]If true, enables printing of the OQ tag with the original base qualities (with -BQSR)")
        #FixMisencodedQuals= define.boolean( name = "Fix Misencoded Quals", default = False, description = "[-fixMisencodedQuals]Fix mis-encoded base quality scores")
        FixMisencodedQuals= define.enum( name = "Fix Misencoded Quals", default = 'False', values = [( 'True','True','' ),( 'False','False','' ),('Auto','Auto','')],  description = "[-fixMisencodedQuals]Fix mis-encoded base quality scores, select Auto to auto detect and fix")
        IntervalMerging= define.enum( name = "Interval Merging", default = "ALL", values = [( 'ALL','ALL','' ),( 'OVERLAPPING_ONLY','OVERLAPPING_ONLY','' )], description = "[-im]Indicates the interval merging rule we should use for abutting intervals")
        IntervalPadding= define.integer( name = "Interval Padding", default = 0, description = "[-ip]Indicates how many basepairs of padding to include around each of the intervals specified with the -L/--intervals argument")
        IntervalSetRule= define.enum( name = "Interval Set Rule", default = "UNION", values = [( 'UNION','UNION','' ),( 'INTERSECTION','INTERSECTION','' )], description = "[-isr]Indicates the set merging approach the interval parser should use to combine the various -L or -XL inputs")
        KeepProgramRecords= define.boolean( name = "Keep Program Records", default = False, description = "[-kpr]Should we override the Walker's default and keep program records from the SAM header")
        MaxRuntime= define.integer( name = "Max Runtime", default = -1, description = "[-maxRuntime]If provided, that GATK will stop execution cleanly as soon after maxRuntime has been exceeded, truncating the run but not exiting with a failure.  By default the value is interpreted in minutes, but this can be changed by maxRuntimeUnits")
        MaxRuntimeUnits= define.enum( name = "Max Runtime Units", default = "MINUTES", values = [( 'NANOSECONDS','NANOSECONDS','' ),( 'MICROSECONDS','MICROSECONDS','' ),( 'MILLISECONDS','MILLISECONDS','' ),( 'SECONDS','SECONDS','' ),( 'MINUTES','MINUTES','' ),( 'HOURS','HOURS','' ),( 'DAYS','DAYS','' )], description = "[-maxRuntimeUnits] The TimeUnit for maxRuntime")
        NonDeterministicRandomSeed= define.boolean( name = "Non Deterministic Random Seed", default = False, description = "[-ndrs]Makes the GATK behave non deterministically, that is, the random numbers generated will be different in every run")
        PedigreeString= define.string( name = "Pedigree String", description = "[-pedString]Pedigree string for samples")
        PedigreeValidationType= define.enum( name = "Pedigree Validation Type", default = "STRICT", values = [( 'STRICT','STRICT','' ),( 'SILENT','SILENT','' )], description = "[-pedValidationType]How strict should we be in validating the pedigree information?")
        PhoneHome= define.enum( name = "Phone Home", default = "STANDARD", values = [( 'NO_ET','NO_ET','' ),( 'STANDARD','STANDARD','' )], description = "[-et]What kind of GATK run report should we generate? STANDARD is the default, can be NO_ET so nothing is posted to the run repository. Please see http://gatkforums.broadinstitute.org/discussion/1250/what-is-phone-home-and-how-does-it-affect-me#latest for details.")
        PreserveQscoresLessThan= define.integer( name = "Preserve Qscores Less Than", default = 6, description = "[-preserveQ]Bases with quality scores less than this threshold won't be recalibrated (with -BQSR)")
        ReadFilter= define.string( name = "Read Filter",default='BadCigar', description = "[-rf]Specify filtration criteria to apply to each read individually")
        ReadGroupBlackList= define.string( name = "Read Group Black List", description = "[-rgbl]Filters out read groups matching : or a .txt file containing the filter strings one per line.")
        RemoveProgramRecords= define.boolean( name = "Remove Program Records", default = False, description = "[-rpr]Should we override the Walker's default and remove program records from the SAM header")
        Tag= define.string( name = "Tag", description = "[-tag]Arbitrary tag string to identify this GATK run as part of a group of runs, for later analysis")
        Unsafe= define.enum( name = "Unsafe", default = "null", values = [( 'ALLOW_UNINDEXED_BAM','ALLOW_UNINDEXED_BAM','' ),( 'ALLOW_UNSET_BAM_SORT_ORDER','ALLOW_UNSET_BAM_SORT_ORDER','' ),( 'NO_READ_ORDER_VERIFICATION','NO_READ_ORDER_VERIFICATION','' ),( 'ALLOW_SEQ_DICT_INCOMPATIBILITY','ALLOW_SEQ_DICT_INCOMPATIBILITY','' ),( 'LENIENT_VCF_PROCESSING','LENIENT_VCF_PROCESSING','' ),( 'ALL','ALL','' ),('null','null','')], description = "[-U]If set, enables unsafe operations: nothing will be checked at runtime.  For expert users only who know what they are doing.  We do not support usage of this argument.")
        UseLegacyDownsampler= define.boolean( name = "Use Legacy Downsampler", default = False, description = "Use the legacy downsampling implementation instead of the newer, less-tested implementation")
        UseOriginalQualities= define.boolean( name = "Use Original Qualities", default = False, description = "[-OQ]If set, use the original base quality scores from the OQ tag when present instead of the standard scores")
        ValidationStrictness= define.enum( name = "Validation Strictness", default = "SILENT", values = [( 'SILENT','SILENT','' ),( 'LENIENT','LENIENT','' ),( 'STRICT','STRICT','' )], description = "[-S]How strict should we be with validation")
        #Groupby= define.enum( name = "Group by", default = "sample", values = [( 'sample_group','sample_group','' ),( 'sample','sample','' ),( 'library','library','' ),( 'platform_unit','platform_unit','' ),( 'chunk','chunk','' ),( 'interval','interval','' )], description = "Inputs will be grouped by selected value from this category. One output will be generated for each group.")
        DivideByIntervals = define.boolean( name = "Divide By Intervals",default = False, description = "Divide the result by Genomic Intervals" )
        #Memoryperjob= define.integer( name = "Memory per job", default = 0, description = "Amount of RAM memory to be used per job. Defaults to 2048MB for Single threaded jobs,and all of the available memory on the instance for multi-threaded jobs. Set to 0 for the default value")
        #Threadsperjob= define.integer( name = "Threads per job", default = 0, description = "For tools which support multiprocessing, this value can be used to set the number of threads to be used. Set to 0 for auto-detect (use with caution,as auto-detect will find the optimal value in most cases)")
        Maximumintervalsize= define.integer( name = "Maximum interval size", default = 500, description = "[-maxInterval]Maximum interval size. Because the realignment algorithm is N^2, allowing too large an interval might take too long to completely realign.")
        Minimumreadsatlocus= define.integer( name = "Minimum reads at locus", default = 4, description = "[-minReads]Minimum reads at a locus to enable using the entropy calculation.")
        Mismatchfraction= define.real( name = "Mismatch fraction", default = 0, description = "[-mismatch]Fraction of base qualities needing to mismatch for a position to have high entropy. To disable this behavior, set this value to <= 0 or > 1. This feature is really only necessary when using an ungapped aligner (e.g. MAQ in the case of single-end read data) and should be used in conjunction with USE_SW' option.")
        Windowsize= define.integer( name = "Window size", default = 10, description = "[-window]Window size for calculating entropy or SNP clusters. Any two SNP calls and/or high entropy positions are considered clustered when they occur no more than this many base pairs apart.")

    def execute(self):
        options=[]
        if(self.params.DisableRandomization):
        #boolean#[-ndrs]Completely eliminates randomization from nondeterministic methods. To be used mostly in the testing framework where dynamic parallelism can result in differing numbers of calls to the generator.
            options.append( '-ndrs')
        if(self.params.AllowPotentiallyMisencodedQuals):
        #boolean#[-allowPotentiallyMisencodedQuals] Do not fail when encountered base qualities that are too high and seemingly indicate a problem with the base quality encoding of the BAM file.
            options.append( '-allowPotentiallyMisencodedQuals')
        if(self.params.BAQCalculationType):
        #enum#[-baq]Type of BAQ calculation to apply in the engine.
            options.extend([ '-baq', self.params.BAQCalculationType])
        if(self.params.BAQGapOpenPenalty):
        #float#[-baqGOP]BAQ gap open penalty (Phred Scaled). Default value is 40. 30 is perhaps better for whole genome call sets
            options.extend([ '-baqGOP', str(self.params.BAQGapOpenPenalty)])
        if(self.params.DefaultBaseQualities):
        #integer#If reads are missing some or all base quality scores, this value will be used for all base quality scores
            options.extend([ '-DBQ',str(self.params.DefaultBaseQualities)])
        if(self.params.DisableIndelQuals):
        #boolean#[-DBQ]If 'true', disables printing of base insertion and base deletion tags (with -BQSR). Turns off printing of the base insertion and base deletion tags when using the -BQSR argument and only the base substitution qualities will be produced.
            options.append( '-DIQ')
        if(self.params.DownsampletoCoverage):
        #integer#[-dcov]Coverage to downsample to at any given locus; note that downsampled reads are randomly selected from all possible reads at a locus. For non-locus-based traversals (eg., ReadWalkers), this sets the maximum number of reads at each alignment start position.
            options.extend([ '-dcov', self.params.DownsampletoCoverage])
        if(self.params.DownsampletoFraction):
        #float#[-dfrac]Fraction [0.0-1.0] of reads to downsample to
            options.extend([ '-dfrac', str(self.params.DownsampletoFraction)])
        if(self.params.DownsamplingType and self.params.DownsamplingType!='null'):
        #enum#[-dt]Type of reads downsampling to employ at a given locus. Reads will be selected randomly to be removed from the pile based on the method described here
            options.extend([ '-dt', self.params.DownsamplingType])
        if(self.params.EmitOriginalQuals):
        #boolean#[-EOQ]If true, enables printing of the OQ tag with the original base qualities (with -BQSR)
            options.append( '-EOQ')
        if(self.params.FixMisencodedQuals == 'True'):
        #boolean#[-fixMisencodedQuals]Fix mis-encoded base quality scores
            options.append( '-fixMisencodedQuals')
        elif ( self.params.FixMisencodedQuals == 'Auto' ):
            if "_quality_scale" in self.inputs.bam[0].meta:
                if ( self.inputs.bam[0].meta.get('_quality_scale') =='Phred+64' ):
                    options.append( '-fixMisencodedQuals')
                    for x in self.inputs.bam:
                        x.meta['_quality_scale'] = 'Phred+33' 
            
        if(self.params.IntervalMerging):
        #enum#[-im]Indicates the interval merging rule we should use for abutting intervals
            options.extend([ '-im', self.params.IntervalMerging])
        if(self.params.IntervalPadding):
        #integer#[-ip]Indicates how many basepairs of padding to include around each of the intervals specified with the -L/--intervals argument
            options.extend([ '-ip', self.params.IntervalPadding])
        if(self.params.IntervalSetRule):
        #enum#[-isr]Indicates the set merging approach the interval parser should use to combine the various -L or -XL inputs
            options.extend([ '-isr', self.params.IntervalSetRule])
        if(self.params.KeepProgramRecords):
        #boolean#[-kpr]Should we override the Walker's default and keep program records from the SAM header
            options.append( '-kpr')
        if(self.params.MaxRuntime):
        #integer#[-maxRuntime]If provided, that GATK will stop execution cleanly as soon after maxRuntime has been exceeded, truncating the run but not exiting with a failure.  By default the value is interpreted in minutes, but this can be changed by maxRuntimeUnits
            options.extend([ '-maxRuntime', self.params.MaxRuntime])
        if(self.params.MaxRuntimeUnits):
        #enum#[-maxRuntimeUnits] The TimeUnit for maxRuntime
            options.extend([ '-maxRuntimeUnits', self.params.MaxRuntimeUnits])
        if(self.params.NonDeterministicRandomSeed):
        #boolean#[-ndrs]Makes the GATK behave non deterministically, that is, the random numbers generated will be different in every run
            options.append( '-ndrs')
        if(self.params.PedigreeString):
        #string#[-pedString]Pedigree string for samples
            options.extend([ '-pedString', self.params.PedigreeString])
        if(self.params.PedigreeValidationType):
        #enum#[-pedValidationType]How strict should we be in validating the pedigree information?
            options.extend([ '-pedValidationType', self.params.PedigreeValidationType])
        if(self.params.PhoneHome):
        #enum#[-et]What kind of GATK run report should we generate? STANDARD is the default, can be NO_ET so nothing is posted to the run repository. Please see http://gatkforums.broadinstitute.org/discussion/1250/what-is-phone-home-and-how-does-it-affect-me#latest for details.
            options.extend([ '-et', self.params.PhoneHome])
        if(self.params.PreserveQscoresLessThan):
        #integer#[-preserveQ]Bases with quality scores less than this threshold won't be recalibrated (with -BQSR)
            options.extend([ '-preserveQ', self.params.PreserveQscoresLessThan])
        if(self.params.ReadFilter):
        #string#[-rf]Specify filtration criteria to apply to each read individually
            options.extend([ '-rf', self.params.ReadFilter])
        if(self.params.ReadGroupBlackList):
        #string#[-rgbl]Filters out read groups matching : or a .txt file containing the filter strings one per line.
            options.extend([ '-rgbl', self.params.ReadGroupBlackList])
        if(self.params.RemoveProgramRecords):
        #boolean#[-rpr]Should we override the Walker's default and remove program records from the SAM header
            options.append( '-rpr')
        if(self.params.Tag):
        #string#[-tag]Arbitrary tag string to identify this GATK run as part of a group of runs, for later analysis
            options.extend([ '-tag', self.params.Tag])
        if(self.params.Unsafe and self.params.Unsafe !='null'):
        #enum#[-U]If set, enables unsafe operations: nothing will be checked at runtime.  For expert users only who know what they are doing.  We do not support usage of this argument.
            options.extend([ '-U', self.params.Unsafe])
        if(self.params.UseLegacyDownsampler):
            options.extend(['-use_legacy_downsampler',self.params.UseLegacyDownsampler])
        #boolean#Use the legacy downsampling implementation instead of the newer, less-tested implementation
        if(self.params.UseOriginalQualities):
        #boolean#[-OQ]If set, use the original base quality scores from the OQ tag when present instead of the standard scores
            options.append( '-OQ')
        if(self.params.ValidationStrictness):
        #enum#[-S]How strict should we be with validation
            options.extend([ '-S', self.params.ValidationStrictness])
#        if(self.params.Groupby):
        #enum#Inputs will be grouped by selected value from this category. One output will be generated for each group.
        #if(self.params.Memoryperjob):
        #integer#Amount of RAM memory to be used per job. Defaults to 2048MB for Single threaded jobs,and all of the available memory on the instance for multi-threaded jobs. Set to 0 for the default value
        #if(self.params.Threadsperjob):
        #integer#For tools which support multiprocessing, this value can be used to set the number of threads to be used. Set to 0 for auto-detect (use with caution,as auto-detect will find the optimal value in most cases)
        if(self.params.Maximumintervalsize):
        #integer#[-maxInterval]Maximum interval size. Because the realignment algorithm is N^2, allowing too large an interval might take too long to completely realign.
            options.extend([ '-maxInterval', self.params.Maximumintervalsize])
        if(self.params.Minimumreadsatlocus):
        #integer#[-minReads]Minimum reads at a locus to enable using the entropy calculation.
            options.extend([ '-minReads', self.params.Minimumreadsatlocus])
        if(self.params.Mismatchfraction):
        #float#[-mismatch]Fraction of base qualities needing to mismatch for a position to have high entropy. To disable this behavior, set this value to <= 0 or > 1. This feature is really only necessary when using an ungapped aligner (e.g. MAQ in the case of single-end read data) and should be used in conjunction with USE_SW' option.
            options.extend([ '-mismatch', str(self.params.Mismatchfraction)])
        if(self.params.Windowsize):
        #integer#[-window]Window size for calculating entropy or SNP clusters. Any two SNP calls and/or high entropy positions are considered clustered when they occur no more than this many base pairs apart.
            options.extend([ '-window', self.params.Windowsize])
        
        


        out_file_name = "dedup.bam.intervals"
        fileNamePath, fileExtension = os.path.splitext( self.inputs.bam[0] )
        out_file_name = fileNamePath + ".dedup.bam.intervals"

    #===========================================================================
    #         # build bam list file
    #     bam_list_file = "bam.list"
    #     with open(bam_list_file, 'w') as f:
    #         for i in range(len(self.inputs.bam_list)/2):
    #             os.rename(self.inputs.bam_list[i*2], self.inputs.bam_list[i*2] + ".bam")
    #             os.rename(self.inputs.bam_list[i*2+1], self.inputs.bam_list[i*2] + ".bai")
    #             f.write("%s\n" % (self.inputs.bam_list[i*2] + ".bam"))
    # 
    #===========================================================================

        if (self.inputs.Gatk_key):
            options.extend(['-K',self.inputs.Gatk_key])
        if (self.inputs.exclude_intervals):
            for x in self.inputs.exclude_intervals:
                options.extend(['-XL',x])
        
        if (self.inputs.exome_bed and not self.params.DivideByIntervals ):
            for x in self.inputs.exome_bed:
                options.extend(['-L',x])
        #run_cmd = ['java', '-Xmx56g', '-Djava.io.tmpdir=/extra/tmp', '-jar', '/opt/bin/GenomeAnalysisTK.jar', '-T', 'RealignerTargetCreator','-nt','32']

        run_cmd = [  '-Djava.io.tmpdir=/extra/tmp', '-jar', '/opt/bin/GenomeAnalysisTK.jar', '-T', 'RealignerTargetCreator'  ]
        run_cmd.extend(options)
        
        #Process('samtools','index',self.inputs.bam).run()
        
        #for  x in self.inputs.bam:
        #    Process('samtools','index',x).run()
        #    run_cmd.extend(['-I',x])
        if  (  self.inputs.bai ):
            run_touch = ['touch']
            for x in self.inputs.bai :
                run_touch.append(x)
            Process(*run_touch).run()
            
        with open('somefile_temp_Samtool_Index.txt', 'a') as the_file:
            for  x in self.inputs.bam:
                fileNamePath2, fileExtension = os.path.splitext(x )
                if ( not ( os.path.exists(fileNamePath2+'.bai') or  os.path.exists(fileNamePath2+'.bam.bai') ) ):
                    the_file.write('samtools index '+ x +'\n')
                run_cmd.extend(['-I',x])
        Process('/opt/bin/multi_process', "-c", '25', '-i', "somefile_temp_Samtool_Index.txt").run()
        
        run_cmd.extend([ '-R', '/opt/db/human_g1k_v37_decoy.fasta',  '--known', '/opt/db/Mills_and_1000G_gold_standard.indels.b37.sites.vcf', '--known', '/opt/db/1000G_phase1.indels.b37.vcf'])
        counter_i = 0
        if (self.params.DivideByIntervals  ):
            with open('somefile_temp_RealignerTargetCreator.txt', 'a') as the_file:
                content=[]
            #    for  Exome_x in self.inputs.exome_bed :
                with open('/opt/db/human_g1k_v37_decoy.breakpoints.bed') as f:
                    content.extend( f.readlines() )
                


                for line in content:
                    run_cmd2 = ['java', '-Xmx2g']
                    run_cmd2.extend(run_cmd)
                    tempstr =  re.split('\t',line)[0]  + ':' + re.split('\t',line)[1] + '-' + re.split('\t',line)[2] 
                    if (tempstr[-1]=='\n'):
                        tempstr = tempstr[:-1]
                    
                    #run_cmd2.extend(['-nt','4','-o', fileNamePath  + ".dedup.bam" + str(counter_i)+ ".intervals",'-L', tempstr ])
                    run_cmd2.extend([ '-o', fileNamePath  + ".dedup.bam." + str(counter_i)+ ".intervals",'-L', tempstr ])

                    the_file.write(' '.join(str(x) for x in run_cmd2) )
                    the_file.write( "\n" )
                    
                    self.outputs.dedup_bam_intervals.add_file ( fileNamePath + ".dedup.bam." + str(counter_i)+ ".intervals")
                    counter_i = counter_i + 1
                    self.outputs.dedup_bam_intervals[-1].meta = self.inputs.bam[0].make_metadata(_interval=tempstr)
                        
            #self.outputs.dedup_bam_intervals[0].meta = self.inputs.bam[0].make_metadata()
                        
            Process("/opt/bin/multi_process",'-c','30' ,"-i",'somefile_temp_RealignerTargetCreator.txt').run()
        
        else:
            run_cmd.extend(  [ '-o', out_file_name ,'-nt','8' ] )
            run_cmd2 = ['java', '-Xmx16g']
            run_cmd2.extend( run_cmd )
            Process(*run_cmd2).run()
            #Process('java', '-Xmx24g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK-2.1-9.jar', '-T', 'RealignerTargetCreator', '-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', '-nt', '12', '-R', '/opt/db/ucsc.hg19.fasta', '-I', bam_list_file, '-o', out_file_name, '--known', '/opt/db/Mills_and_1000G_gold_standard.indels.hg19.vcf', '--known', '/opt/db/1000G_phase1.indels.hg19.vcf', '-rf', 'BadCigar', '-L', self.inputs.exome_bed).run()
        #else:
        #    Process('java', '-Xmx24g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK-2.1-9.jar', '-T', 'RealignerTargetCreator', '-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', '-nt', '12', '-R', '/opt/db/ucsc.hg19.fasta', '-I', bam_list_file, '-o', out_file_name, '--known', '/opt/db/Mills_and_1000G_gold_standard.indels.hg19.vcf', '--known', '/opt/db/1000G_phase1.indels.hg19.vcf', '-rf', 'BadCigar').run()
            self.outputs.dedup_bam_intervals.add_file ( out_file_name )
            self.outputs.dedup_bam_intervals[-1].meta = self.inputs.bam[0].make_metadata()
        
def test_wrapperclass():
        #inputs = {'bam': ['/l3bioinfo/test-data/ERR315384_1.fastq.clip.bam'],'exome_bed':'/l3bioinfo/test-data/human_g1k_v37_decoy.breakpoints.bed' }
        inputs = {'bam': ['/l3bioinfo/test-data/ERR174325_1.fastq.bam'],'bai':['/l3bioinfo/test-data/ERR174325_1.fastq.bai'],'exome_bed':['/l3bioinfo/test-data/exome_targets.b37.bed'] }
        #inputs = {'bam': '/l3bioinfo/test-data/ERR17432_1.fastq.clip.bam','exome_bed':'/l3bioinfo/test-data/exome_targets.b37.bed' }
        #inputs = {'bam': ['/l3bioinfo/test-data/ERR174324_1.fastq.clip.bam','/l3bioinfo/test-data/ERR174325_1.fastq.clip.bam','/l3bioinfo/test-data/ERR174326_1.fastq.clip.bam','/l3bioinfo/test-data/ERR174327_1.fastq.clip.bam'],'bai':  ['/l3bioinfo/test-data/ERR174324_1.fastq.clip.bam.bai','/l3bioinfo/test-data/ERR174325_1.fastq.clip.bam.bai','/l3bioinfo/test-data/ERR174326_1.fastq.clip.bam.bai','/l3bioinfo/test-data/ERR174327_1.fastq.clip.bam.bai'],'exome_bed':'/l3bioinfo/test-data/human_g1k_v37_decoy.breakpoints.bed' }
        params = {}
        outputs = WrapperClass(inputs, params).test()
        #assert outputs.out.endswith('mock.bam')
