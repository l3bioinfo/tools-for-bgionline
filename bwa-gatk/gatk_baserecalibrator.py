import os
from l3sdk import define, Process, require

@require(mem_mb=18000,high_io=True,disk_space_gb=150)
class WrapperClass(define.Wrapper):
    class Inputs(define.Inputs):
        bam = define.input(name="Read sequence", required=True,description = "Read sequence in BAM format", list=True )  #, alt_path='/extra'
        bai = define.input(name="BAM index files", required=False,description="According BAM index",list=True )
        exclude_intervals = define.input(name="Exclude Intervals",list=True, required=False)
        exome_bed = define.input(name="Genomic Intervals",list=True, required=False)
        Gatk_key =  define.input(name="GATK key", required=False)

    class Outputs(define.Outputs):
        #dedup_realn_bam_grp = define.output(description="Dedup Realn BAM Group", required=True)
        grp = define.output(name="Recalibration Table Output", description="The output recalibration table file to create", required=True, list=False  )  #, alt_path='/extra'
    class Params(define.Params):
        BinaryTagName= define.string( name = "Binary Tag Name", description = "[-bintag] the binary tag covariate name if using it")
        Covariate= define.enum( name = "Covariate", values = [], description = "[-cov]One or more covariates to be used in the recalibration. Can be specified multiple times")
        DeletionsDefaultQuality= define.integer( name = "Deletions Default Quality", default = 45, description = "[-ddq] default quality for the base deletions covariate")
        IndelsContextSize= define.integer( name = "Indels Context Size", default = 3, description = "[-ics] size of the k-mer context to be used for base insertions and deletions")
        Disableindelquals= define.boolean( name = "Disable indel quals", default = True, description = "[--disable_indel_quals] Disable indel quality recalibration. Must be set to true in GATK Lite.")
        InsertionsDefaultQuality= define.integer( name = "Insertions Default Quality", default = 45, description = "[-idq] default quality for the base insertions covariate")
        LowQualityTail= define.integer( name = "Low Quality Tail", default = 2, description = "[-lqt] minimum quality for the bases in the tail of the reads to be considered")
        MaximumCycleValue= define.integer( name = "Maximum Cycle Value", default = 500, description = "[-maxCycle ] the maximum cycle value permitted for the Cycle covariate")
        MismatchesContextSize= define.integer( name = "Mismatches Context Size", default = 2, description = "[-maxCycle]size of the k-mer context to be used for base mismatches")
        MismatchesDefaultQuality= define.integer( name = "Mismatches Default Quality", default = -1, description = "[-mdq]default quality for the base mismatches covariate")
        NoStandardCovs= define.boolean( name = "No Standard Covs", default = False, description = "[-noStandard]Do not use the standard set of covariates, but rather just the ones listed using the -cov argument. Cannot be used if grouped by interval.")
        QuantizingLevels= define.integer( name = "Quantizing Levels", default = 16, description = "[-ql] number of distinct quality scores in the quantized output")
        SolidNocallStrategy= define.enum( name = "Solid Nocall Strategy", default = "THROW_EXCEPTION", values = [( 'THROW_EXCEPTION','THROW_EXCEPTION','' ),( 'LEAVE_READ_UNRECALIBRATED','LEAVE_READ_UNRECALIBRATED','' ),( 'PURGE_READ','PURGE_READ','' )], description = "[--solid_nocall_strategy]Defines the behavior of the recalibrator when it encounters no calls in the color space. Options = THROW_EXCEPTION, LEAVE_READ_UNRECALIBRATED, or PURGE_READ")
        SolidRecalMode= define.enum( name = "Solid Recal Mode", default = "SET_Q_ZERO", values = [( 'DO_NOTHING','DO_NOTHING','' ),( 'SET_Q_ZERO','SET_Q_ZERO','' ),( 'SET_Q_ZERO_BASE_N','SET_Q_ZERO_BASE_N','' ),( 'REMOVE_REF_BIAS','REMOVE_REF_BIAS','' )], description = "[-sMode]How should we recalibrate solid bases in which the reference was inserted? Options = DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, or REMOVE_REF_BIAS")
        BqsrBaqGapOpenPenalty= define.real( name = "Bqsr Baq Gap Open Penalty", default = 40, description = "[-bqsrBAQGOP]BQSR BAQ gap open penalty (Phred Scaled).  Default value is 40.  30 is perhaps better for whole genome call sets")
        RunWithoutDbsnpPotentiallyRuiningQuality= define.boolean( name = "Run Without Dbsnp Potentially Ruining Quality", default = False, description = "[-run_without_dbsnp_potentially_ruining_quality] If specified, allows the recalibrator to be used without a dbsnp rod. Very unsafe and for expert users only.")


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
        FixMisencodedQuals= define.boolean( name = "Fix Misencoded Quals", default = False, description = "[-fixMisencodedQuals]Fix mis-encoded base quality scores")
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
        #Memoryperjob= define.integer( name = "Memory per job", default = 0, description = "Amount of RAM memory to be used per job. Defaults to 2048MB for Single threaded jobs,and all of the available memory on the instance for multi-threaded jobs. Set to 0 for the default value")
        #Threadsperjob= define.integer( name = "Threads per job", default = 0, description = "For tools which support multiprocessing, this value can be used to set the number of threads to be used. Set to 0 for auto-detect (use with caution,as auto-detect will find the optimal value in most cases)")

    def execute(self):
        options=[]

        if (self.params.BinaryTagName):
        #string[-bintag] the binary tag covariate name if using it"
            options.extend(['-bintag', self.params.BinaryTagName])
        if (self.params.Covariate):
        #enum[-cov]One or more covariates to be used in the recalibration. Can be specified multiple times"
            options.extend(['-cov', self.params.Covariate])
        if (self.params.DeletionsDefaultQuality):
        #integer[-ddq] default quality for the base deletions covariate
            options.extend(['-ddq', self.params.DeletionsDefaultQuality])
        if (self.params.IndelsContextSize):
        #integer[-ics] size of the k-mer context to be used for base insertions and deletions
            options.extend(['-ics', self.params.IndelsContextSize])
        if (self.params.InsertionsDefaultQuality):
        #integer[-idq] default quality for the base insertions covariate
            options.extend(['-idq', self.params.InsertionsDefaultQuality])
        if (self.params.LowQualityTail):
        #integer[-lqt] minimum quality for the bases in the tail of the reads to be considered
            options.extend(['-lqt', self.params.LowQualityTail])
        if (self.params.MaximumCycleValue):
        #integer[-maxCycle ] the maximum cycle value permitted for the Cycle covariate
            options.extend(['-maxCycle', self.params.MaximumCycleValue])
        if (self.params.MismatchesContextSize):
        #integer= "[-mcs]size of the k-mer context to be used for base mismatches")
            options.extend(['-mcs', self.params.MismatchesContextSize])
        if (self.params.MismatchesDefaultQuality):
        #integer[-mdq]default quality for the base mismatches covariate
            options.extend(['-mdq', self.params.MismatchesDefaultQuality])
        if (self.params.QuantizingLevels):
        #integer[-ql] number of distinct quality scores in the quantized output
            options.extend(['-ql', self.params.QuantizingLevels])
        if (self.params.SolidNocallStrategy):
        #enum[--solid_nocall_strategy]Defines the behavior of the recalibrator when it encounters no calls in the color space. Options = THROW_EXCEPTION, LEAVE_READ_UNRECALIBRATED, or PURGE_READ
            options.extend(['--solid_nocall_strategy', self.params.SolidNocallStrategy])
        if (self.params.SolidRecalMode):
        #enum[-sMode]How should we recalibrate solid bases in which the reference was inserted? Options = DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, or REMOVE_REF_BIAS
            options.extend(['-sMode', self.params.SolidRecalMode])
        if (self.params.BqsrBaqGapOpenPenalty):
        #real[-bqsrBAQGOP]BQSR BAQ gap open penalty (Phred Scaled).  Default value is 40.  30 is perhaps better for whole genome call sets
            options.extend(['-bqsrBAQGOP', self.params.BqsrBaqGapOpenPenalty])
        if (self.params.RunWithoutDbsnpPotentiallyRuiningQuality):
        #boolean[-run_without_dbsnp_potentially_ruining_quality] If specified, allows the recalibrator to be used without a dbsnp rod. Very unsafe and for expert users only.
            options.append('-run_without_dbsnp_potentially_ruining_quality')

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
        if(self.params.FixMisencodedQuals):
        #boolean#[-fixMisencodedQuals]Fix mis-encoded base quality scores
            options.append( '-fixMisencodedQuals')
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

        


        out_file_name = "dedup.realn.bam.table"
        fileNamePath, fileExtension = os.path.splitext( self.inputs.bam[0] )
        out_file_name = fileNamePath + ".dedup.realn.bam.table"
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
        if  (  self.inputs.bai ):
            run_touch = ['touch']
            for x in self.inputs.bai :
                run_touch.append(x)
            Process(*run_touch).run()
        
        if (self.inputs.Gatk_key):
            options.extend(['-K',self.inputs.Gatk_key])
        if (self.inputs.exclude_intervals):
            for x in self.inputs.exclude_intervals:
                options.extend(['-XL',x])
                # 16g
        #run_cmd = ['java', '-Xmx56g', '-Djava.io.tmpdir=/extra/tmp', '-jar', '/opt/bin/GenomeAnalysisTK.jar', '-T', 'BaseRecalibrator','-nct', '32' ,'--disable_indel_quals']
        run_cmd = ['java', '-Xmx16g', '-Djava.io.tmpdir=/extra/tmp', '-jar', '/opt/bin/GenomeAnalysisTK.jar', '-T', 'BaseRecalibrator','-nct', '8' ,'--disable_indel_quals']
        run_cmd.extend(options)
        with open('somefile_temp_Samtool_Index.txt', 'a') as the_file:
            for i in xrange(0,len(self.inputs.bam)):
                fileNamePath2, fileExtension = os.path.splitext( self.inputs.bam[i] )
                if ( not ( os.path.exists(fileNamePath2+'.bai') or  os.path.exists(fileNamePath2+'.bam.bai') ) ):
                    the_file.write('samtools index '+ self.inputs.bam[i] +'\n')
        Process('/opt/bin/multi_process', "-c", '25', '-i', "somefile_temp_Samtool_Index.txt").run()
        run_cmd.extend([ '-R', '/opt/db/human_g1k_v37_decoy.fasta', '-o', out_file_name ])
        for i in xrange (0,len(self.inputs.bam)):
            run_cmd.extend(['-I', self.inputs.bam[i]])
        run_cmd.extend(['--knownSites','/opt/db/dbsnp_137.b37.vcf' ])
        #run_cmd.extend(['--intervals', '/opt/db/BaseRecali.intervals']);
        for x in self.inputs.exome_bed:
            run_cmd.extend(['-L',x])
                
        #Process('echo', *run_cmd).run()
        Process(*run_cmd).run()

        

    #===========================================================================
    # if self.params.fix_misencoded_quality_scores:
    #     Process('java', '-Xmx16g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK-3.2-2.jar', '-et', 'NO_ET', '-K', '/opt/db/rbluo_cs.hku.hk.key', '-T', 'BaseRecalibrator', '--knownSites', '/opt/db/Mills_and_1000G_gold_standard.indels.hg19.vcf', '--knownSites', '/opt/db/1000G_phase1.indels.hg19.vcf', '-nct', '12', '--knownSites', '/opt/db/dbsnp_138.hg19.vcf', '-R', '/opt/db/ucsc.hg19.fasta', '-I', in0, '-I', in1, '-I', in2, '-I', in3, '-I', in4, '-I', in5, '-o', out_file_name, '-rf', 'BadCigar', '-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', '-L', 'chr1', '-L', 'chr6', '-L', 'chr22', '-L', 'chrX', '--fix_misencoded_quality_scores').run()
    # else:
    #     Process('java', '-Xmx16g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK-3.2-2.jar', '-et', 'NO_ET', '-K', '/opt/db/rbluo_cs.hku.hk.key', '-T', 'BaseRecalibrator', '--knownSites', '/opt/db/Mills_and_1000G_gold_standard.indels.hg19.vcf', '--knownSites', '/opt/db/1000G_phase1.indels.hg19.vcf', '-nct', '12', '--knownSites', '/opt/db/dbsnp_138.hg19.vcf', '-R', '/opt/db/ucsc.hg19.fasta', '-I', in0, '-I', in1, '-I', in2, '-I', in3, '-I', in4, '-I', in5, '-o', out_file_name, '-rf', 'BadCigar', '-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', '-L', 'chr1', '-L', 'chr6', '-L', 'chr22', '-L', 'chrX').run()
    #===========================================================================

        self.outputs.grp = out_file_name
        self.outputs.grp.meta = self.inputs.bam.make_metadata()
        
def test_wrapperclass():
        #inputs = {'bam': '/l3bioinfo/test-data/ERR315384_1.fastq.clip.dedup.realn.bam','exome_bed':'/l3bioinfo/test-data/exome_targets.b37.bed' }
        inputs = {'bam': { '/l3bioinfo/test-data/ERR315384_1.fastq.clip.dedup.realn.bam', '/l3bioinfo/test-data/ERR315384_2.fastq.clip.dedup.realn.bam'},'exome_bed':'/l3bioinfo/test-data/exome_targets.b37.bed' }
        #inputs = {'bam': '/l3bioinfo/test-data/L3_HUM.dedup.realn.bam','exome_bed':'/l3bioinfo/test-data/human_g1k_v37_decoy.breakpoints.bed' }
        params = {}
        #params = {'AllowPotentiallyMisencodedQuals':True,'FixMisencodedQuals':True}
        outputs = WrapperClass(inputs, params).test()
