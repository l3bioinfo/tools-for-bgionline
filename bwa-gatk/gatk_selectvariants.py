import os, re
from l3sdk import define, Process, require

@require(mem_mb=3000)
class WrapperClass(define.Wrapper):
    class Inputs(define.Inputs):
        inp = define.input(description="Input VCF", required=True)

        exclude_intervals = define.input(name="Exclude Intervals",list=True, required=False)
        exome_bed = define.input(name="Genomic Intervals",list=True, required=False)
        Gatk_key =  define.input(name="GATK key", required=False)
        
        concordance = define.input(name="Concordance", required=False)
        discordance = define.input(name="Discordance", required=False)
        keepIDs = define.input(name="Keep IDs", required=False)
    class Outputs(define.Outputs):
        out = define.output(description="SNP Output VCF" )

    class Params(define.Params):
        AllowNonoverlappingCommandLineSamples= define.boolean( name = "Allow Nonoverlapping Command Line Samples", default = False, description = "[--ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES]Allow a samples other than those in the VCF to be specified on the command line. These samples will be ignored.")
        ExcludeSampleName= define.string( name = "Exclude Sample Name", description = "[--exclude_sample_name] Exclude genotypes from this sample. Can be specified multiple times")
        ExcludeFiltered= define.boolean( name = "Exclude Filtered", default = False, description = "[-ef]Don't include filtered loci in the analysis")
        ExcludeNonVariants= define.boolean( name = "Exclude Non Variants", default = False, description = "[-env]Don't include loci found to be non-variant after the subsetting procedure")
        KeepOriginalAc= define.boolean( name = "Keep Original Ac", default = False, description = "[--keepOriginalAC]Store the original AC, AF, and AN values in the INFO field after selecting (using keys AC_Orig, AF_Orig, and AN_Orig)")
        MaxIndelSize= define.integer( name = "Max Indel Size", default = 2147483647, description = "[--maxIndelSize]indel size select")
        MendelianViolation= define.boolean( name = "Mendelian Violation", default = False, description = "[-mv]output mendelian violation sites only")
        Mvq= define.real( name = "Mvq", default = 0, description = "[-mvq]Minimum genotype QUAL score for each trio member required to accept a site as a violation")
        Regenotype= define.boolean( name = "Regenotype", default = False, description = "[-regenotype]re-genotype the selected samples based on their GLs (or PLs)")
        RemoveFractionGenotypes= define.real( name = "Remove Fraction Genotypes", default = 0, description = "[-fractionGenotypes]Selects a fraction (a number between 0 and 1) of the total genotypes at random from the variant track and sets them to nocall")
        RestrictAllelesTo= define.enum( name = "Restrict Alleles To", default = "ALL", values = [( 'ALL','ALL','' ),( 'MULTIALLELIC','MULTIALLELIC','' ),( 'BIALLELIC','BIALLELIC','' )], description = "[--restrictAllelesTo]Select only variants of a particular allelicity. Valid options are ALL (default), MULTIALLELIC or BIALLELIC")
        SampleExpressions= define.string( name = "Sample Expressions", description = "[-se]Regular expression to select many samples from the ROD tracks provided. Can be specified multiple times")
        SampleName= define.string( name = "Sample Name", description = "[-sn]Include genotypes from this sample. Can be specified multiple times")
        SelectExpressions= define.string( name = "Select Expressions", description = "[-select]One or more criteria to use when selecting the data")
        SelectRandomFraction= define.real( name = "Select Random Fraction", default = 0, description = "[-fraction]Selects a fraction (a number between 0 and 1) of the total variants at random from the variant track")
        SelectTypeToInclude= define.string( name = "Select Type To Include",  description = "[-selectType] Select only a certain type of variants from the input file. Valid types are INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION. List them out and separate by comma or space")
        
        
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

        
        if(self.params.AllowNonoverlappingCommandLineSamples):
        #boolean#[--ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES]Allow a samples other than those in the VCF to be specified on the command line. These samples will be ignored.
            options.append( '--ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES')
        if(self.params.ExcludeSampleName):
        #string#[--exclude_sample_name] Exclude genotypes from this sample. Can be specified multiple times
            options.extend([ '--exclude_sample_name', self.params.ExcludeSampleName])
        if(self.params.ExcludeFiltered):
        #boolean#[-ef]Don't include filtered loci in the analysis
            options.append( '-ef')
        if(self.params.ExcludeNonVariants):
        #boolean#[-env]Don't include loci found to be non-variant after the subsetting procedure
            options.append( '-env')
        if(self.params.KeepOriginalAc):
        #boolean#[--keepOriginalAC]Store the original AC, AF, and AN values in the INFO field after selecting (using keys AC_Orig, AF_Orig, and AN_Orig)
            options.append( '--keepOriginalAC')
        if(self.params.MaxIndelSize):
        #integer#[--maxIndelSize]indel size select
            options.extend([ '--maxIndelSize', str(self.params.MaxIndelSize)])
        if(self.params.MendelianViolation):
        #boolean#[-mv]output mendelian violation sites only
            options.append( '-mv')
        if(self.params.Mvq):
        #float#[-mvq]Minimum genotype QUAL score for each trio member required to accept a site as a violation
            options.extend([ '-mvq', str(self.params.Mvq)])
        if(self.params.Regenotype):
        #boolean#[-regenotype]re-genotype the selected samples based on their GLs (or PLs)
            options.append( '-regenotype')
        if(self.params.RemoveFractionGenotypes):
        #float#[-fractionGenotypes]Selects a fraction (a number between 0 and 1) of the total genotypes at random from the variant track and sets them to nocall
            options.extend([ '-fractionGenotypes', str(self.params.RemoveFractionGenotypes)])
        if(self.params.RestrictAllelesTo):
        #enum#[--restrictAllelesTo]Select only variants of a particular allelicity. Valid options are ALL (default), MULTIALLELIC or BIALLELIC
            options.extend([ '--restrictAllelesTo', self.params.RestrictAllelesTo])
        if(self.params.SampleExpressions):
        #string#[-se]Regular expression to select many samples from the ROD tracks provided. Can be specified multiple times
            options.extend([ '-se', self.params.SampleExpressions])
        if(self.params.SampleName):
        #string#[-sn]Include genotypes from this sample. Can be specified multiple times
            options.extend([ '-sn', self.params.SampleName])
        if(self.params.SelectExpressions):
        #string#[-select]One or more criteria to use when selecting the data
            options.extend([ '-select', self.params.SelectExpressions])
        if(self.params.SelectRandomFraction):
        #float#[-fraction]Selects a fraction (a number between 0 and 1) of the total variants at random from the variant track
            options.extend([ '-fraction', str(self.params.SelectRandomFraction)])

        
        #if(self.params.SelectTypeToInclude):
        #enum#[-selectType] Select only a certain type of variants from the input file. Valid types are INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION. Can be specified multiple times
        #    options.extend([ '-selectType', self.params.SelectTypeToInclude])
        if (self.params.SelectTypeToInclude):
            for x in re.split(',| ' , self.params.SelectTypeToInclude):
                options.extend(['-selectType',str(x)])  

        out_file_name = "snp_out.vcf"
        fileNamePath, fileExtension = os.path.splitext( self.inputs.inp )
        out_file_name = fileNamePath + ".snp_out.vcf"
        #Process('java', '-Xmx2g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK-3.2-2.jar', '-R', '/opt/db/ucsc.hg19.fasta', '-et', 'NO_ET', '-K', '/opt/db/rbluo_cs.hku.hk.key', '-T', 'SelectVariants', '--variant', self.inputs.inp, '-selectType', 'SNP', '-o', out_file_name).run()
        if (self.inputs.Gatk_key):
            options.extend(['-K',self.inputs.Gatk_key])
        if (self.inputs.exclude_intervals):
            for x in self.inputs.exclude_intervals:
                options.extend(['-XL',x])
        if (self.inputs.exome_bed):
            for x in self.inputs.exome_bed:
                options.extend(['-L',x])

        if (self.inputs.concordance):
            options.extend(['--concordance',self.inputs.concordance])                
        if (self.inputs.discordance):
            options.extend(['--discordance',self.inputs.discordance])                
        if (self.inputs.keepIDs):
            options.extend(['--keepIDs',self.inputs.keepIDs])
                            
        run_cmd = ['java', '-Xmx2g', '-Djava.io.tmpdir=/tmp', '-jar', '/opt/bin/GenomeAnalysisTK.jar', '-T', 'SelectVariants']
        run_cmd.extend(options)
        
        run_cmd.extend([ '-R', '/opt/db/human_g1k_v37_decoy.fasta', '--variant', self.inputs.inp, '-o', out_file_name])
        Process(*run_cmd).run()    
        
        
        
        
        self.outputs.out = out_file_name
        self.outputs.out.meta = self.inputs.inp.make_metadata()

def test_wrapperclass():
        inputs = {'inp': '/l3bioinfo/test-data/all.vcf' }
        params = { 'SelectExpressions':'QD < 2.0 && MD < 40.0 && FS > 60.0 && HaplotypeScore > 13.0 && MQRankSum < -12.5 && ReadPosRankSum < -8.0'}
        outputs = WrapperClass(inputs, params).test()

