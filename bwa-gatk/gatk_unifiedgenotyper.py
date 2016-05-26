import os
from l3sdk import define, Process, require

@require(mem_mb=18000, high_io=True,disk_space_gb=150)
class WrapperClass(define.Wrapper):
    class Inputs(define.Inputs):
        bam = define.input(name="Read sequence", required=True, list=True, description = "Read sequence in BAM format"   )
        bai = define.input(name="BAM index files", required=False,description="According BAM index",list=True)
        exclude_intervals = define.input(name="Exclude Intervals",list=True, required=False)
        exome_bed = define.input(name="Genomic Intervals",list=True, required=False)
        Gatk_key =  define.input(name="GATK key", required=False)

        dbSNP  =  define.input(name="dbSNP", required=False)
        Alleles  =  define.input(name="Alleles", required=False)
        comp  =  define.input(name="comp", required=False, description = 'Comparison file in VCF format')
        #ref_sample_calls = define.input(name="Reference Sample Calls", required=False, description = 'VCF file with the truth callset for the reference sample')
        BQSR = define.input(name="BQSR", required=False, description = 'Input covariates table file for on-the-fly base quality score recalibration') 
    class Outputs(define.Outputs):
        all_vcf = define.output(description="Output VCF", list=True, required=True )  #, alt_path='/extra'

    class Params(define.Params):
        Annotation= define.string( name = "Annotation", description = "[-A]One or more specific annotations to apply to variant calls")
        ComputeSlod= define.boolean( name = "Compute Slod", default = False, description = "[-slod]If provided, we will calculate the SLOD (SB annotation)")
        Contamination= define.real( name = "Contamination", default = 0.05, description = "[-contamination]Fraction of contamination in sequencing data (for all samples) to aggressively remove.")
        ExcludeAnnotation= define.string( name = "Exclude Annotation", description = "[-XA]One or more specific annotations to exclude")
        GenotypeLikelihoodsModel= define.enum( name = "Genotype Likelihoods Model", default = "SNP", values = [( 'SNP','SNP','' ),( 'INDEL','INDEL','' ),( 'GENERALPLOIDYSNP','GENERALPLOIDYSNP','' ),( 'GENERALPLOIDYINDEL','GENERALPLOIDYINDEL','' ),( 'BOTH','BOTH','' )], description = "[-glm]Genotype likelihoods calculation model to employ -- SNP is the default option, while INDEL is also available for calling indels and BOTH is available for calling both together")
        GenotypingMode= define.enum( name = "Genotyping Mode", default = "DISCOVERY", values = [( 'DISCOVERY','DISCOVERY','' ),( 'GENOTYPE_GIVEN_ALLELES','GENOTYPE_GIVEN_ALLELES','' )], description = "[-gt_mode]Specifies how to determine the alternate alleles to use for genotyping")
        Group= define.string( name = "Group", default = "Standard", description = "[-G]One or more classes/groups of annotations to apply to variant calls")
        Heterozygosity= define.real( name = "Heterozygosity", default = 0.001, description = "[-hets] Heterozygosity value used to compute prior likelihoods for any locus")
        IgnoreLaneInfo= define.boolean( name = "Ignore Lane Info", default = False, description = "[-ignoreLane] Ignore lane when building error model, error model is then per-site")
        IndelHeterozygosity= define.real( name = "Indel Heterozygosity", default = 0.000125, description = "[-indelHeterozygosity]Heterozygosity for indel calling")
        MaxDeletionFraction= define.real( name = "Max Deletion Fraction", default = 0.05, description = "[-deletions]Maximum fraction of reads with deletions spanning this locus for it to be callable [to disable, set to  1; default:0.05]")
        MinBaseQualityScore= define.integer( name = "Min Base Quality Score", default = 17, description = "[-mbq]Minimum base quality required to consider a base for calling")
        MinIndelCnt= define.integer( name = "Min Indel Cnt", default = 5, description = "[-minIndelCnt]Minimum number of consensus indels required to trigger genotyping run")
        MinIndelFrac= define.real( name = "Min Indel Frac", default = 0.25, description = "[-minIndelFrac]Minimum fraction of all reads at a locus that must contain an indel (of any allele) for that sample to contribute to the indel count for alleles")
        OutputMode= define.enum( name = "Output Mode", default = "EMIT_VARIANTS_ONLY", values = [( 'EMIT_VARIANTS_ONLY','EMIT_VARIANTS_ONLY','' ),( 'EMIT_ALL_CONFIDENT_SITES','EMIT_ALL_CONFIDENT_SITES','' ),( 'EMIT_ALL_SITES','EMIT_ALL_SITES','' )], description = "[-out_mode]Specifies which type of calls we should output")
        PairHmmImplementation= define.enum( name = "Pair Hmm Implementation", default = "ORIGINAL", values = [( 'EXACT','EXACT','' ),( 'ORIGINAL','ORIGINAL','' ),( 'CACHING','CACHING','' ),( 'LOGLESS_CACHING','LOGLESS_CACHING','' )], description = "[-pairHMM]The PairHMM implementation to use for -glm INDEL genotype likelihood calculations")
        PcrErrorRate= define.real( name = "Pcr Error Rate", default = 0.0001, description = "The PCR error rate to be used for computing fragment-based likelihoods")
        StandCallConf= define.real( name = "Stand Call Conf", default = 30, description = "[-stand_call_conf]The minimum phred-scaled confidence threshold at which variants should be called")
        StandEmitConf= define.real( name = "Stand Emit Conf", default = 30, description = "[-stand_emit_conf]The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold)")
        IndelGapContinuationPenalty= define.integer( name = "Indel Gap Continuation Penalty", default = 10, description = "[-indelGCP]Indel gap continuation penalty, as Phred-scaled probability.  I.e., 30 => 10^-30/10")
        IndelGapOpenPenalty= define.integer( name = "Indel Gap Open Penalty", default = 45, description = "[-indelGOP]Indel gap open penalty, as Phred-scaled probability.  I.e., 30 => 10^-30/10")
        MaxAlternateAlleles= define.integer( name = "Max Alternate Alleles", default = 6, description = "[-maxAltAlleles]Maximum number of alternate alleles to genotype")
        PNonrefModel= define.enum( name = "P Nonref Model", default = "EXACT_INDEPENDENT", values = [( 'EXACT_INDEPENDENT','EXACT_INDEPENDENT','' ),( 'EXACT_REFERENCE','EXACT_REFERENCE','' ),( 'EXACT_ORIGINAL','EXACT_ORIGINAL','' ),( 'EXACT_GENERAL_PLOIDY','EXACT_GENERAL_PLOIDY','' )], description = "[--p_nonref_model] Non-reference probability calculation model to employ")
        
        
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
        DivideByIntervals = define.boolean ( name = "Divide By Intervals, if true, it assumes '_interval' is set for each input bam file ", default = False, description = "Divide the result by Genomic Intervals" )
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
     
            
    
        if(self.params.Annotation):
        #string#[-A]One or more specific annotations to apply to variant calls
            options.extend([ '-A', self.params.Annotation])
        if(self.params.ComputeSlod):
        #boolean#[-slod]If provided, we will calculate the SLOD (SB annotation)
            options.append( '-slod')
        if(self.params.Contamination):
        #float#[-contamination]Fraction of contamination in sequencing data (for all samples) to aggressively remove.
            options.extend([ '-contamination', str(self.params.Contamination)])
        if(self.params.ExcludeAnnotation):
        #string#[-XA]One or more specific annotations to exclude
            options.extend([ '-XA', self.params.ExcludeAnnotation])
        if(self.params.GenotypeLikelihoodsModel):
        #enum#[-glm]Genotype likelihoods calculation model to employ -- SNP is the default option, while INDEL is also available for calling indels and BOTH is available for calling both together
            options.extend([ '-glm', self.params.GenotypeLikelihoodsModel])
        if(self.params.GenotypingMode):
        #enum#[-gt_mode]Specifies how to determine the alternate alleles to use for genotyping
            options.extend([ '-gt_mode', self.params.GenotypingMode])
        if(self.params.Group):
        #string#[-G]One or more classes/groups of annotations to apply to variant calls
            options.extend([ '-G', self.params.Group])
        if(self.params.Heterozygosity):
        #float#[-hets] Heterozygosity value used to compute prior likelihoods for any locus
            options.extend([ '-hets', str(self.params.Heterozygosity)])
        if(self.params.IgnoreLaneInfo):
        #boolean#[-ignoreLane] Ignore lane when building error model, error model is then per-site
            options.append( '-ignoreLane')
        if(self.params.IndelHeterozygosity):
        #float#[-indelHeterozygosity]Heterozygosity for indel calling
            options.extend([ '-indelHeterozygosity', str(self.params.IndelHeterozygosity)])
        if(self.params.MaxDeletionFraction):
        #float#[-deletions]Maximum fraction of reads with deletions spanning this locus for it to be callable [to disable, set to  1; default:0.05]
            options.extend([ '-deletions', str(self.params.MaxDeletionFraction)])
        if(self.params.MinBaseQualityScore):
        #integer#[-mbq]Minimum base quality required to consider a base for calling
            options.extend([ '-mbq', str(self.params.MinBaseQualityScore)])
        if(self.params.MinIndelCnt):
        #integer#[-minIndelCnt]Minimum number of consensus indels required to trigger genotyping run
            options.extend([ '-minIndelCnt', str(self.params.MinIndelCnt)])
        if(self.params.MinIndelFrac):
        #float#[-minIndelFrac]Minimum fraction of all reads at a locus that must contain an indel (of any allele) for that sample to contribute to the indel count for alleles
            options.extend([ '-minIndelFrac', str(self.params.MinIndelFrac)])
        if(self.params.OutputMode):
        #enum#[-out_mode]Specifies which type of calls we should output
            options.extend([ '-out_mode', self.params.OutputMode])
        if(self.params.PairHmmImplementation):
        #enum#[-pairHMM]The PairHMM implementation to use for -glm INDEL genotype likelihood calculations
            options.extend([ '-pairHMM', self.params.PairHmmImplementation])
        if(self.params.PcrErrorRate):
        #float#The PCR error rate to be used for computing fragment-based likelihoods
            options.extend(['--pcr_error_rate', str(self.params.PcrErrorRate) ])
        if(self.params.StandCallConf):
        #float#[-stand_call_conf]The minimum phred-scaled confidence threshold at which variants should be called
            options.extend([ '-stand_call_conf', str(self.params.StandCallConf)])
        if(self.params.StandEmitConf):
        #float#[-stand_emit_conf]The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold)
            options.extend([ '-stand_emit_conf', str(self.params.StandEmitConf)])
        if(self.params.IndelGapContinuationPenalty):
        #integer#[-indelGCP]Indel gap continuation penalty, as Phred-scaled probability.  I.e., 30 => 10^-30/10
            options.extend([ '-indelGCP', str(self.params.IndelGapContinuationPenalty)])
        if(self.params.IndelGapOpenPenalty):
        #integer#[-indelGOP]Indel gap open penalty, as Phred-scaled probability.  I.e., 30 => 10^-30/10
            options.extend([ '-indelGOP', str(self.params.IndelGapOpenPenalty)])
        if(self.params.MaxAlternateAlleles):
        #integer#[-maxAltAlleles]Maximum number of alternate alleles to genotype
            options.extend([ '-maxAltAlleles', str(self.params.MaxAlternateAlleles)])
        if(self.params.PNonrefModel):
        #enum#[--pnrm] Non-reference probability calculation model to employ
            options.extend([ '--p_nonref_model', self.params.PNonrefModel])
            
        out_file_name = "all.vcf"
        fileNamePath, fileExtension = os.path.splitext( self.inputs.bam[0] )
        out_file_name = fileNamePath + ".vcf"
        #out_file_name = os.path.basename(self.inputs.bam)
        #out_file_name = outFileNameVCF[::-1].replace(".bam"[::-1], ".vcf"[::-1], 1)[::-1]

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
                
        if (self.inputs.Alleles):
            options.extend(['--alleles',self.inputs.Alleles])
        if (self.inputs.comp):
            options.extend(['--comp',self.inputs.comp])
        if (self.inputs.dbSNP):
            options.extend(['--dbsnp',self.inputs.dbSNP])
        if (self.inputs.BQSR):
            options.extend(['--BQSR',self.inputs.BQSR])
            
        if  (  self.inputs.bai ):
            run_touch = ['touch']
            for x in self.inputs.bai :
                run_touch.append(x)
            Process(*run_touch).run() 
            
        with open('somefile_temp_Samtool_Index.txt', 'a') as the_file:
            for i in xrange(0,len(self.inputs.bam)):
                fileNamePath2, fileExtension = os.path.splitext( self.inputs.bam[i] )
                if ( not ( os.path.exists(fileNamePath2+'.bai') or  os.path.exists(fileNamePath2+'.bam.bai') ) ):
                    the_file.write('samtools index '+ self.inputs.bam[i] +'\n')
        
        Process('/opt/bin/multi_process', "-c", '25', '-i', "somefile_temp_Samtool_Index.txt").run()
        
        if ( self.params.DivideByIntervals ):
            with open('somefile_temp_UnifiedGenotyper.txt', 'a') as the_file:
                for i in xrange(0,len(self.inputs.bam)):
                    fileNamePath, fileExtension = os.path.splitext( self.inputs.bam[i] )
                    
                    run_cmd2 = ['java','-Xmx2g', '-Djava.io.tmpdir=/extra/tmp', '-jar', '/opt/bin/GenomeAnalysisTK.jar', '-T', 'UnifiedGenotyper', '-nt', '4', '-nct', '1']
                    run_cmd2.extend(options)
                    run_cmd2.extend(['--dbsnp', '/opt/db/dbsnp_137.b37.vcf'])
                    run_cmd2.extend(['-R','/opt/db/human_g1k_v37_decoy.fasta'])
                    run_cmd2.extend(['-I', self.inputs.bam[i], '-L', self.inputs.bam[i].meta.get('_interval')])
                    run_cmd2.extend(['-o', '%s.vcf' % (fileNamePath, )])
                    the_file.write(' '.join( str(x) for x in run_cmd2 ) + '\n' )

            Process('/opt/bin/multi_process', "-c", '25', '-i', "somefile_temp_UnifiedGenotyper.txt").run()
            for i in xrange(0,len(self.inputs.bam)):
                fileNamePath, fileExtension = os.path.splitext( self.inputs.bam[i] )
                
                self.outputs.all_vcf.add_file( '%s.vcf' % ( fileNamePath,  ) )
                self.outputs.all_vcf[-1].meta = self.inputs.bam[i].make_metadata()
                # remove _interval meta
                if '_interval' in self.inputs.bam[i].meta:
                    self.outputs.all_vcf[-1].meta.pop("_interval", None)
        else:
            #run_cmd = ['java', '-Xmx56g', '-Djava.io.tmpdir=/extra/tmp', '-jar', '/opt/bin/GenomeAnalysisTK.jar', '-T', 'UnifiedGenotyper','-nt', '32', '-nct', '1']
            run_cmd = ['java', '-Xmx16g', '-Djava.io.tmpdir=/extra/tmp', '-jar', '/opt/bin/GenomeAnalysisTK.jar', '-T', 'UnifiedGenotyper','-nt', '12', '-nct', '1']
            run_cmd.extend(options)
            run_cmd.extend(['--dbsnp','/opt/db/dbsnp_137.b37.vcf' ])
            run_cmd.extend([ '-R', '/opt/db/human_g1k_v37_decoy.fasta', '-I', self.inputs.bam[0], '-o', out_file_name])
            #Process('samtools', 'index', self.inputs.bam[0]).run()
            Process(*run_cmd).run()
            self.outputs.all_vcf.add_file(out_file_name)
            self.outputs.all_vcf[-1].meta = self.inputs.bam[0].make_metadata()

def test_wrapperclass():
        inputs = {'bam': {'/l3bioinfo/test-data/output1.bam', '/l3bioinfo/test-data/output2.bam' } }
        params = {'DivideByIntervals': True }
        outputs = WrapperClass(inputs, params).test()
