import os
import re
from l3sdk import define, Process, require

@require(mem_mb=5000,high_io=True,disk_space_gb=150)
class WrapperClass(define.Wrapper):
    class Inputs(define.Inputs):
        bam = define.input(description="BAM File", required=True ,  list=True )
        bai = define.input(name="BAM index files", required=False,description="According BAM index",list=True )
        dedup_bam_intervals = define.input(description="Dedup BAM Intervals", list=True, required=True)
        exclude_intervals = define.input(name="Exclude Intervals",list=True, required=False)
        exome_bed = define.input(name="Genomic Intervals",list=True, required=False)
        Gatk_key =  define.input(name="GATK key", required=False)
        
    class Outputs(define.Outputs):
        dedup_realn_bam = define.output(description="Dedup Realn BAM", required=True, list=True  )
        dedup_realn_bai = define.output(description="Dedup Realn BAM Index", required=True, list=True  )

    class Params(define.Params):
        ConsensusDeterminationModel= define.enum( name = "Consensus Determination Model", default = "USE_READS", values = [( 'KNOWNS_ONLY','KNOWNS_ONLY','' ),( 'USE_READS','USE_READS','' ),( 'USE_SW','USE_SW','' )], description = "[-model]Determines how to compute the possible alternate consenses")
        LodThresholdForCleaning= define.real( name = "Lod Threshold For Cleaning", default = 5, description = "[-LOD]LOD threshold above which the cleaner will clean")
        EntropyThreshold= define.real( name = "Entropy Threshold", default = 0.15, description = "[-entropy]percentage of mismatches at a locus to be considered having high entropy")
        MaxConsensuses= define.integer( name = "Max Consensuses", default = 30, description = "[--maxConsensuses]max alternate consensuses to try (necessary to improve performance in deep coverage)")
        MaxIsizeForMovement= define.integer( name = "Max Isize For Movement", default = 3000, description = "[-maxIsize] maximum insert size of read pairs that we attempt to realign")
        MaxPositionalMoveAllowed= define.integer( name = "Max Positional Move Allowed", default = 200, description = "[-maxPosMove]maximum positional move in basepairs that a read can be adjusted during realignment")
        MaxReadsForConsensuses= define.integer( name = "Max Reads For Consensuses", default = 120, description = "[-greedy] max reads used for finding the alternate consensuses (necessary to improve performance in deep coverage)")
        MaxReadsForRealignment= define.integer( name = "Max Reads For Realignment", default = 20000, description = "[-maxReads]max reads allowed at an interval for realignment")
        MaxReadsInMemory= define.integer( name = "Max Reads In Memory", default = 150000, description = "[-maxInMemory]max reads allowed to be kept in memory at a time by the SAMFileWriter")
        NoOriginalAlignmentTags= define.boolean( name = "No Original Alignment Tags", default = False, description = "[-noTags]Don't output the original cigar or alignment start tags for each realigned read in the output bam")

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

        fileNamePath, fileExtension = os.path.splitext( self.inputs.bam[0] )
        out_file_name = fileNamePath + ".indel.bam"

        #out_file_prefix = "indel"
        #out_file_name = out_file_prefix + ".bam"
    # build bam list file
    
        if (self.inputs.Gatk_key):
            options.extend(['-K',self.inputs.Gatk_key])
        if (self.inputs.exclude_intervals):
            for x in self.inputs.exclude_intervals:
                options.extend(['-XL',x])
        if (self.inputs.exome_bed):
            for x in self.inputs.exome_bed:
                options.extend(['-L',x])
                
        #run_cmd = ['java', '-Xmx4g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK.jar', '-T', 'IndelRealigner' ]

        run_cmd = [ '-Djava.io.tmpdir=/extra/tmp', '-jar', '/opt/bin/GenomeAnalysisTK.jar', '-T', 'IndelRealigner' ]
        run_cmd.extend(options)
        
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
        
        run_cmd.extend([ '-R', '/opt/db/human_g1k_v37_decoy.fasta' ])
        counter_i = 0
        if (self.params.DivideByIntervals  ):
            with open('somefile_temp_IndelRealigner.txt', 'a') as the_file:
                #content=[]
                #for  Exome_x in self.inputs.exome_bed :
                #    with open(Exome_x) as f:
                #        content.extend( f.readlines() )


                    # for line in content:
                for interval_x in self.inputs.dedup_bam_intervals :
                    # y = [int(x.group()) for x in re.finditer(r'\d+', interval_x ) ]
                    # counter_i = int( y[-1] )
                    # line = content [ counter_i ]

                    run_cmd2 = ['java', '-Xmx2g']
                    run_cmd2.extend(run_cmd)
                    run_cmd2.extend(['--targetIntervals',interval_x])
                    run_cmd2.extend(['-o', fileNamePath + ".indel"+ str(counter_i)+".bam"  ])

                    self.outputs.dedup_realn_bam.add_file ( fileNamePath + ".indel."+ str(counter_i)+".bam" )
                    self.outputs.dedup_realn_bai.add_file ( fileNamePath + ".indel."+ str(counter_i)+".bai" )
                    if '_interval' in interval_x.meta:
                        temp_interval = interval_x.meta.get('_interval')
                        run_cmd2.extend(  [ '-L', temp_interval  ] )
                        self.outputs.dedup_realn_bam[-1].meta = self.inputs.bam[0].make_metadata( _interval=temp_interval )
                    else:
                        self.outputs.dedup_realn_bam[-1].meta = self.inputs.bam[0].make_metadata()
                        
                    the_file.write(' '.join(str(x) for x in run_cmd2) )
                    the_file.write( "\n" )
                    
                    counter_i = counter_i + 1
     
            Process("/opt/bin/multi_process",'-c','30' ,"-i",'somefile_temp_IndelRealigner.txt').run()
        
        else:
            if (self.inputs.dedup_bam_intervals):
                for x in self.inputs.dedup_bam_intervals:
                    run_cmd.extend(['--targetIntervals', x ])

            run_cmd.extend(  [ '-o', out_file_name  ] )
            run_cmd2 = ['java', '-Xmx4g']
            run_cmd2.extend( run_cmd )
            Process(*run_cmd2).run()
            #Process('java', '-Xmx24g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK-2.1-9.jar', '-T', 'RealignerTargetCreator', '-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', '-nt', '12', '-R', '/opt/db/ucsc.hg19.fasta', '-I', bam_list_file, '-o', out_file_name, '--known', '/opt/db/Mills_and_1000G_gold_standard.indels.hg19.vcf', '--known', '/opt/db/1000G_phase1.indels.hg19.vcf', '-rf', 'BadCigar', '-L', self.inputs.exome_bed).run()
        #else:
        #    Process('java', '-Xmx24g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK-2.1-9.jar', '-T', 'RealignerTargetCreator', '-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', '-nt', '12', '-R', '/opt/db/ucsc.hg19.fasta', '-I', bam_list_file, '-o', out_file_name, '--known', '/opt/db/Mills_and_1000G_gold_standard.indels.hg19.vcf', '--known', '/opt/db/1000G_phase1.indels.hg19.vcf', '-rf', 'BadCigar').run()
            self.outputs.dedup_realn_bam.add_file ( out_file_name )
            self.outputs.dedup_realn_bai.add_file ( fileNamePath + ".indel.bai" )
            self.outputs.dedup_realn_bam[-1].meta = self.inputs.bam[0].make_metadata()
        

        
        
        
        
        
#        run_cmd = ['java', '-Xmx56g', '-Djava.io.tmpdir=/extra/tmp', '-jar', '/opt/bin/GenomeAnalysisTK.jar', '-T', 'IndelRealigner' ]
#        run_cmd.extend(options)
#        Process('samtools','index',self.inputs.bam).run()

#        run_cmd.extend([ '-R', '/opt/db/human_g1k_v37_decoy.fasta', '-I', self.inputs.bam, '-o', out_file_name ])
#        Process(*run_cmd).run()

        #chr_list = self.params.chromosomes.split(' ')

       # Process('java', '-Xmx4g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK-2.1-9.jar', '-et', 'NO_ET', '-K', '/opt/db/rbluo_cs.hku.hk.key', '-T', 'IndelRealigner', '-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', '-model', 'USE_READS', '-known', '/opt/db/Mills_and_1000G_gold_standard.indels.hg19.vcf', '-known', '/opt/db/1000G_phase1.indels.hg19.vcf', '-R', '/opt/db/ucsc.hg19.fasta', '--targetIntervals', self.inputs.dedup_bam_intervals, '-I', bam_list_file, '-o', out_file_name, '-rf', 'BadCigar', *chr_list).run()

#        self.outputs.dedup_realn_bam =  out_file_name
#        self.outputs.dedup_realn_bam.meta = self.inputs.bam.make_metadata()
        
def test_wrapperclass():
        temp_list = []
        for x in range(0,100):
            temp_list.append( '/l3bioinfo/test-data/ERR315384_1.fastq.clip.dedup.bam'+str(x)+'.intervals')
        inputs = {'bam': ['/l3bioinfo/test-data/ERR315384_1.fastq.clip.bam'],'dedup_bam_intervals': temp_list ,'exome_bed':'/l3bioinfo/test-data/human_g1k_v37_decoy.breakpoints.bed'  }
        #inputs = {'bam': ['/l3bioinfo/test-data/ERR315384_1.fastq.clip.bam'],'dedup_bam_intervals': ['/l3bioinfo/test-data/ERR315384_1.fastq.clip.dedup.bam.intervals'] ,'exome_bed':'/l3bioinfo/test-data/human_g1k_v37_decoy.breakpoints.bed'  }
        params = {'DivideByIntervals':True}
        outputs = WrapperClass(inputs, params).test()
        #assert outputs.out.endswith('mock.bam')
