import os
import re
from l3sdk import define, Process, require

@require(mem_mb=3000)
class WrapperClass(define.Wrapper):
    class Inputs(define.Inputs):
        vcf = define.input(description="Input VCF", required=True)
        exclude_intervals = define.input(name="Exclude Intervals",list=True, required=False)
        exome_bed = define.input(name="Genomic Intervals",list=True, required=False)
        Gatk_key =  define.input(name="GATK key", required=False)

        mask =  define.input(name="Mask", required=False,description ='Input ROD mask')

    class Outputs(define.Outputs):
        out = define.output(name = 'output VCF file', description="Filtered VCF" )

    class Params(define.Params):
        ClusterSize= define.integer( name = "Cluster Size", default = 3, description = "[-cluster]The number of SNPs which make up a cluster")
        ClusterWindowSize= define.integer( name = "Cluster Window Size", default = 0, description = "[-window]The window size (in bases) in which to evaluate clustered SNPs")
        FiltersName = define.string( name = "FilterNames", description = "[--filterName] Names to use for the list of filters, divide them using comma")
        Filters= define.string( name = "Filters", description = "[-filter] One or more expression used with INFO fields to filter, divide them using comma")
        GenotypefiltersName = define.string( name = "Genotype filters Names", description ="[--genotypeFilterName] Names to use for the list of sample/genotype filters (must be a 1-to-1 mapping); this name is put in the FILTER field for variants that get filtered, divide them using comma" )
        Genotypefilters= define.string( name = "Genotype filters", description = "[-G_filter] One or more expression used with FORMAT (sample/genotype-level) fields to filter (see documentation guide for more info), divide them using comma")

        InvalidatePreviousFilters= define.boolean( name = "Invalidate Previous Filters", default = False, description = "[--invalidatePreviousFilters]Remove previous filters applied to the VCF")
        MaskExtension= define.integer( name = "Mask Extension", default = 0, description = "[-maskExtend]How many bases beyond records from a provided 'mask' rod should variants be filtered")
        MaskName= define.string( name = "Mask Name", default = "Mask", description = "[--maskName]The text to put in the FILTER field if a 'mask' rod is provided and overlaps with a variant call")
        MissingValuesInExpressionsShouldEvaluateAsFailing= define.boolean( name = "Missing Values In Expressions Should Evaluate As Failing", default = False, description = "[--missingValuesInExpressionsShouldEvaluateAsFailing] When evaluating the JEXL expressions, missing values should be considered failing the expression")
        
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

            
        if(self.params.ClusterSize):
        #integer#[-cluster]The number of SNPs which make up a cluster
            options.extend([ '-cluster', str(self.params.ClusterSize)])
        if(self.params.ClusterWindowSize):
        #integer#[-window]The window size (in bases) in which to evaluate clustered SNPs
            options.extend([ '-window', str(self.params.ClusterWindowSize)])
        if(self.params.FiltersName):
        #string#[-filter] One or more expression used with INFO fields to filter
            for x in re.split(',',self.params.FiltersName):
                options.extend([ '--filterName', x])
        if(self.params.Filters):
        #string#[-filter] One or more expression used with INFO fields to filter
            for x in re.split(',',self.params.Filters):
                options.extend([ '-filter', x])
        if(self.params.GenotypefiltersName):
        #string#[-G_filter] One or more expression used with FORMAT (sample/genotype-level) fields to filter (see documentation guide for more info)
            for x in re.split(',',self.params.GenotypefiltersName):
                options.extend([ '--genotypeFilterName', x])
        if(self.params.Genotypefilters):
        #string#[-G_filter] One or more expression used with FORMAT (sample/genotype-level) fields to filter (see documentation guide for more info)
            for x in re.split(',',self.params.Genotypefilters):
                options.extend([ '-G_filter', x ])
        if(self.params.InvalidatePreviousFilters):
        #boolean#[--invalidatePreviousFilters]Remove previous filters applied to the VCF
            options.append( '--invalidatePreviousFilters')
        if(self.params.MaskExtension):
        #integer#[-maskExtend]How many bases beyond records from a provided 'mask' rod should variants be filtered
            options.extend([ '-maskExtend', str(self.params.MaskExtension)])
        if(self.params.MaskName):
        #string#[--maskName]The text to put in the FILTER field if a 'mask' rod is provided and overlaps with a variant call
            options.extend([ '--maskName', self.params.MaskName])
        if(self.params.MissingValuesInExpressionsShouldEvaluateAsFailing):
        #boolean#[--missingValuesInExpressionsShouldEvaluateAsFailing] When evaluating the JEXL expressions, missing values should be considered failing the expression
            options.append( '--missingValuesInExpressionsShouldEvaluateAsFailing')


        out_file_name = 'filtered.vcf'
        fileNamePath, fileExtension = os.path.splitext( self.inputs.vcf )
        out_file_name = fileNamePath + ".filtered.vcf"
    #===========================================================================
    # if self.inputs.exome_bed:
    #     Process('java', '-Xmx16g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK-3.2-2.jar', '-l', 'INFO', '-R', '/opt/db/ucsc.hg19.fasta', '-et', 'NO_ET', '-K', '/opt/db/rbluo_cs.hku.hk.key', '-T', 'VariantRecalibrator', '-nt', '12', '-input', self.inputs.inp, '-resource:hapmap,known=false,training=true,truth=true,prior=15.0', '/opt/db/hapmap_3.3.hg19.vcf', '-resource:omni,known=false,training=true,truth=false,prior=12.0', '/opt/db/1000G_omni2.5.hg19.vcf', '-resource:dbsnp,known=true,training=false,truth=false,prior=8.0', '/opt/db/dbsnp_138.hg19.vcf', '-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0', '/opt/db/Mills_and_1000G_gold_standard.indels.hg19.vcf', '-resource:phase1,VCF,known=true,training=true,truth=true,prior=9.0', '/opt/db/1000G_phase1.indels.hg19.vcf', '-an', 'FS', '-an', 'QD', '-an', 'ReadPosRankSum', '-an', 'HaplotypeScore', '-an', 'MQ', '-recalFile', recal_file_name, '-tranchesFile', tranches_file_name, '--TStranche', '90.0', '--TStranche', '93.0', '--TStranche', '95.0', '--TStranche', '97.0', '--TStranche', '99.0', '--TStranche', '100.0', '-mode', 'BOTH', '-L', self.inputs.exome_bed).run()
    # else:
    #     Process('java', '-Xmx16g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK-3.2-2.jar', '-l', 'INFO', '-R', '/opt/db/ucsc.hg19.fasta', '-et', 'NO_ET', '-K', '/opt/db/rbluo_cs.hku.hk.key', '-T', 'VariantRecalibrator', '-nt', '12', '-input', self.inputs.inp, '-resource:hapmap,known=false,training=true,truth=true,prior=15.0', '/opt/db/hapmap_3.3.hg19.vcf', '-resource:omni,known=false,training=true,truth=false,prior=12.0', '/opt/db/1000G_omni2.5.hg19.vcf', '-resource:dbsnp,known=true,training=false,truth=false,prior=8.0', '/opt/db/dbsnp_138.hg19.vcf', '-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0', '/opt/db/Mills_and_1000G_gold_standard.indels.hg19.vcf', '-resource:phase1,VCF,known=true,training=true,truth=true,prior=9.0', '/opt/db/1000G_phase1.indels.hg19.vcf', '-an', 'FS', '-an', 'QD', '-an', 'ReadPosRankSum', '-an', 'HaplotypeScore', '-an', 'MQ', '-recalFile', recal_file_name, '-tranchesFile', tranches_file_name, '--TStranche', '90.0', '--TStranche', '93.0', '--TStranche', '95.0', '--TStranche', '97.0', '--TStranche', '99.0', '--TStranche', '100.0', '-mode', 'BOTH').run()
    #===========================================================================
        if (self.inputs.Gatk_key):
            options.extend(['-K',self.inputs.Gatk_key])
        if (self.inputs.exclude_intervals):
            for x in self.inputs.exclude_intervals:
                options.extend(['-XL',x])
        if (self.inputs.exome_bed):
            for x in self.inputs.exome_bed:
                options.extend(['-L',x])
        if (self.inputs.mask):
            options.extend(['--mask',self.inputs.mask])
        
        run_cmd = ['java', '-Xmx2g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK.jar', '-T', 'VariantFiltration' ]
        run_cmd.extend(options)
        
        run_cmd.extend([ '-R', '/opt/db/human_g1k_v37_decoy.fasta', '--variant', self.inputs.vcf, '-o', out_file_name])
        Process(*run_cmd).run()
        
        self.outputs.out = out_file_name
        self.outputs.out.meta = self.inputs.vcf.make_metadata()

def test_wrapperclass():
        inputs = {'vcf': '/l3bioinfo/test-data/snp_out.vcf' }
        params = {'Filters': 'QD < 2.0,MD < 40.0,FS > 60.0','FiltersName':'QD,MD,FS'}
        outputs = WrapperClass(inputs, params).test()