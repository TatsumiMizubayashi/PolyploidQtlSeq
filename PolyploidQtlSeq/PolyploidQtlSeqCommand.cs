using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeqCore.Application.Pipeline;
using PolyploidQtlSeqCore.Mapping;
using PolyploidQtlSeqCore.VariantCall;
using System.ComponentModel.DataAnnotations;
using System.Reflection;
using chr = PolyploidQtlSeqCore.QtlAnalysis.Chr;
using d = PolyploidQtlSeqCore.QtlAnalysis.Distribution;
using io = PolyploidQtlSeqCore.QtlAnalysis.IO;
using og = PolyploidQtlSeqCore.QtlAnalysis.OxyGraph;
using op = PolyploidQtlSeqCore.Options;
using qa = PolyploidQtlSeqCore.QtlAnalysis;
using qsf = PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter;
using sw = PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeq
{
    /// <summary>
    /// メインコマンド
    /// </summary>
    [Command(Description = "Run the Polyploid QTL-seq pipeline.")]
    [VersionOptionFromMember("-v|--version", Description = "Show version number.", MemberName = nameof(GetVersion))]
    [Subcommand(
        typeof(QualityControlCommand),
        typeof(QtlAnalysisCommand))]
    internal sealed class PolyploidQtlSeqCommand : CommandBase, IQtlSeqPipelineSettingValue
    {
        /// <summary>
        /// 高次倍数性QTL-Seqコマンドを作成する。
        /// </summary>
        public PolyploidQtlSeqCommand()
        {
            ReferenceSequence = "";
            Parent1Dir = "";
            Parent2Dir = "";
            Bulk1Dir = "";
            Bulk2Dir = "";
            ChrSizeThreshold = chr.ChrSizeThreshold.DEFAULT;
            AnalysisChrNames = "";
            MinMq = MinmumMappingQuality.DEFAULT;
            MinBq = MinmumBaseQuality.DEFAULT;
            AdjustMq = AdjustMappingQuality.DEFAULT;
            SnpEffMaxHeap = PolyploidQtlSeqCore.VariantCall.SnpEffMaxHeap.DEFAULT;
            SnpEffConfigFile = "";
            SnpEffDatabaseName = "";
            OutputDir = "";
            Parent1MostAlleleRateThreshold = qsf.Parent1MostAlleleRateThreshold.DEFAULT;
            Parent2SnpIndexRange = qsf.Parent2SnpIndexRange.DEFAULT;
            MinimumDepthThreshold = qsf.MinimumDepthThreshold.DEFAULT;
            MaxBulkSnpIndexThreshold = qsf.MaxBulkSnpIndexThreshold.DEFAULT;
            Ploidy = d.Ploidy.DEFAULT;
            Parent2PlexNumber = d.Parent2PlexNumber.DEFAULT;
            Bulk1Number = d.Bulk1Number.DEFAULT;
            Bulk2Number = d.Bulk2Number.DEFAULT;
            ReplicatesNumber = d.ReplicatesNumber.DEFAULT;
            WindowSize = sw.WindowSize.DEFAULT;
            StepSize = sw.StepSize.DEFAULT;
            FigureWidth = og.FigureWidth.DEFAULT;
            FigureHeight = og.FigureHeight.DEFAULT;
            XAxisMajorStep = og.XAxisMajorStep.DEFAULT;
            DisplayAnnotationImpacts = io.DisplayAnnotationImpacts.DEFAULT;
            ThreadNumber = PolyploidQtlSeqCore.Share.ThreadNumber.DEFAULT;
            ParameterFile = "";
        }

        [Option(ShortName = PolyploidQtlSeqCore.Share.ReferenceSequence.SHORT_NAME, LongName = PolyploidQtlSeqCore.Share.ReferenceSequence.LONG_NAME,
            Description = PolyploidQtlSeqCore.Share.ReferenceSequence.DESCRIPTION, ValueName = "")]
        public string ReferenceSequence { get; set; }

        [Option(ShortName = Parent1Directory.SHORT_NAME, LongName = Parent1Directory.LONG_NAME,
            Description = Parent1Directory.DESCRIPTION, ValueName = "")]
        public string Parent1Dir { get; set; }

        [Option(ShortName = Parent2Directory.SHORT_NAME, LongName = Parent2Directory.LONG_NAME,
            Description = Parent2Directory.DESCRIPTION, ValueName = "")]
        public string Parent2Dir { get; set; }

        [Option(ShortName = Bulk1Directory.SHORT_NAME, LongName = Bulk1Directory.LONG_NAME,
            Description = Bulk1Directory.DESCRIPTION, ValueName = "")]
        public string Bulk1Dir { get; set; }

        [Option(ShortName = Bulk2Directory.SHORT_NAME, LongName = Bulk2Directory.LONG_NAME,
            Description = Bulk2Directory.DESCRIPTION, ValueName = "")]
        public string Bulk2Dir { get; set; }

        [Option(ShortName = qa.OutputDirectory.SHORT_NAME, LongName = qa.OutputDirectory.LONG_NAME,
            Description = qa.OutputDirectory.DESCRIPTION, ValueName = "")]
        public string OutputDir { get; set; }

        [Option(ShortName = chr.ChrSizeThreshold.SHORT_NAME, LongName = chr.ChrSizeThreshold.LONG_NAME,
            Description = chr.ChrSizeThreshold.DESCRIPTION, ValueName = "")]
        [Range(chr.ChrSizeThreshold.MINIMUM, chr.ChrSizeThreshold.MAXIMUM, ErrorMessage = chr.ChrSizeThreshold.VALIDATION_ERROR_MESSAGE)]
        public int ChrSizeThreshold { get; set; }

        [Option(ShortName = chr.AnalysisChrNames.SHORT_NAME, LongName = chr.AnalysisChrNames.LONG_NAME,
            Description = chr.AnalysisChrNames.DESCRIPTION, ValueName = "")]
        public string AnalysisChrNames { get; set; }

        [Option(ShortName = MinmumMappingQuality.SHORT_NAME, LongName = MinmumMappingQuality.LONG_NAME,
            Description = MinmumMappingQuality.DESCRIPTION, ValueName = "")]
        [Range(MinmumMappingQuality.MINIMUM, MinmumMappingQuality.MAXIMUM, ErrorMessage = MinmumMappingQuality.VALIDATION_ERROR_MESSAGE)]
        public int MinMq { get; set; }

        [Option(ShortName = MinmumBaseQuality.SHORT_NAME, LongName = MinmumBaseQuality.LONG_NAME,
            Description = MinmumBaseQuality.DESCRIPTION, ValueName = "")]
        [Range(MinmumBaseQuality.MINIMUM, MinmumBaseQuality.MAXIMUM, ErrorMessage = MinmumBaseQuality.VALIDATION_ERROR_MESSAGE)]
        public int MinBq { get; set; }

        [Option(ShortName = AdjustMappingQuality.SHORT_NAME, LongName = AdjustMappingQuality.LONG_NAME,
            Description = AdjustMappingQuality.DESCRIPTION, ValueName = "")]
        [Range(AdjustMappingQuality.MINIMUM, AdjustMappingQuality.MAXIMUM, ErrorMessage = AdjustMappingQuality.VALIDATION_ERROR_MESSAGE)]
        public int AdjustMq { get; set; }

        [Option(ShortName = PolyploidQtlSeqCore.VariantCall.SnpEffMaxHeap.SHORT_NAME, LongName = PolyploidQtlSeqCore.VariantCall.SnpEffMaxHeap.LONG_NAME,
            Description = PolyploidQtlSeqCore.VariantCall.SnpEffMaxHeap.DESCRIPTION, ValueName = "")]
        [Range(PolyploidQtlSeqCore.VariantCall.SnpEffMaxHeap.MINIMUM, PolyploidQtlSeqCore.VariantCall.SnpEffMaxHeap.MAXIMUM, ErrorMessage = PolyploidQtlSeqCore.VariantCall.SnpEffMaxHeap.VALIDATION_ERROR_MESSAGE)]
        public int SnpEffMaxHeap { get; set; }

        [Option(ShortName = PolyploidQtlSeqCore.VariantCall.SnpEffConfigFile.SHORT_NAME, LongName = PolyploidQtlSeqCore.VariantCall.SnpEffConfigFile.LONG_NAME,
            Description = PolyploidQtlSeqCore.VariantCall.SnpEffConfigFile.DESCRIPTION, ValueName = "")]
        public string SnpEffConfigFile { get; set; }

        [Option(ShortName = SnpEffDatabase.SHORT_NAME, LongName = SnpEffDatabase.LONG_NAME,
            Description = SnpEffDatabase.DESCRIPTION, ValueName = "")]
        public string SnpEffDatabaseName { get; set; }


        [Option(ShortName = qsf.Parent1MostAlleleRateThreshold.SHORT_NAME, LongName = qsf.Parent1MostAlleleRateThreshold.LONG_NAME,
            Description = qsf.Parent1MostAlleleRateThreshold.DESCRIPTION, ValueName = "")]
        public double Parent1MostAlleleRateThreshold { get; set; }

        [Option(ShortName = qsf.Parent2SnpIndexRange.SHORT_NAME, LongName = qsf.Parent2SnpIndexRange.LONG_NAME,
            Description = qsf.Parent2SnpIndexRange.DESCRIPTION, ValueName = "")]
        public string Parent2SnpIndexRange { get; set; }

        [Option(ShortName = qsf.MinimumDepthThreshold.SHORT_NAME, LongName = qsf.MinimumDepthThreshold.LONG_NAME,
            Description = qsf.MinimumDepthThreshold.DESCRIPTION, ValueName = "")]
        public int MinimumDepthThreshold { get; set; }

        [Option(ShortName = qsf.MaxBulkSnpIndexThreshold.SHORT_NAME, LongName = qsf.MaxBulkSnpIndexThreshold.LONG_NAME,
            Description = qsf.MaxBulkSnpIndexThreshold.DESCRIPTION, ValueName = "")]
        public double MaxBulkSnpIndexThreshold { get; set; }


        [Option(ShortName = d.Ploidy.SHORT_NAME, LongName = d.Ploidy.LONG_NAME,
            Description = d.Ploidy.DESCRIPTION, ValueName = "")]
        [Range(d.Ploidy.MINIMUM, d.Ploidy.MAXIMUM, ErrorMessage = d.Ploidy.VALIDATION_ERROR_MESSAGE)]
        public int Ploidy { get; set; }

        [Option(ShortName = d.Parent2PlexNumber.SHORT_NAME, LongName = d.Parent2PlexNumber.LONG_NAME,
            Description = d.Parent2PlexNumber.DESCRIPTION, ValueName = "")]
        [Range(d.Parent2PlexNumber.MINIMUM, d.Parent2PlexNumber.MAXIMUM,
            ErrorMessage = d.Parent2PlexNumber.VALIDATION_ERROR_MESSAGE)]
        public int Parent2PlexNumber { get; set; }

        [Option(ShortName = d.Bulk1Number.SHORT_NAME, LongName = d.Bulk1Number.LONG_NAME,
            Description = d.Bulk1Number.DESCRIPTION, ValueName = "")]
        [Range(d.Bulk1Number.MINIMUM, d.Bulk1Number.MAXIMUM, ErrorMessage = d.Bulk1Number.VALIDATION_ERROR_MESSAGE)]
        public int Bulk1Number { get; set; }

        [Option(ShortName = d.Bulk2Number.SHORT_NAME, LongName = d.Bulk2Number.LONG_NAME,
            Description = d.Bulk2Number.DESCRIPTION, ValueName = "")]
        [Range(d.Bulk2Number.MINIMUM, d.Bulk2Number.MAXIMUM, ErrorMessage = d.Bulk2Number.VALIDATION_ERROR_MESSAGE)]
        public int Bulk2Number { get; set; }

        [Option(ShortName = d.ReplicatesNumber.SHORT_NAME, LongName = d.ReplicatesNumber.LONG_NAME,
            Description = d.ReplicatesNumber.DESCRIPTION, ValueName = "")]
        [Range(d.ReplicatesNumber.MINIMUM, d.ReplicatesNumber.MAXIMUM, ErrorMessage = d.ReplicatesNumber.VALIDATION_ERROR_MESSAGE)]
        public int ReplicatesNumber { get; set; }

        [Option(ShortName = sw.WindowSize.SHORT_NAME, LongName = sw.WindowSize.LONG_NAME,
            Description = sw.WindowSize.DESCRIPTION, ValueName = "")]
        [Range(sw.WindowSize.MINIMUM, sw.WindowSize.MAXIMUM, ErrorMessage = sw.WindowSize.VALIDATION_ERROR_MESSAGE)]
        public int WindowSize { get; set; }

        [Option(ShortName = sw.StepSize.SHORT_NAME, LongName = sw.StepSize.LONG_NAME,
            Description = sw.StepSize.DESCRIPTION, ValueName = "")]
        [Range(sw.StepSize.MINIMUM, sw.StepSize.MAXIMUM, ErrorMessage = sw.StepSize.VALIDATION_ERROR_MESSAGE)]
        public int StepSize { get; set; }

        [Option(ShortName = og.FigureWidth.SHORT_NAME, LongName = og.FigureWidth.LONG_NAME,
            Description = og.FigureWidth.DESCRIPTION, ValueName = "")]
        [Range(og.FigureWidth.MINIMUM, og.FigureWidth.MAXIMUM, ErrorMessage = og.FigureWidth.VALIDATION_ERROR_MESSAGE)]
        public int FigureWidth { get; set; }

        [Option(ShortName = og.FigureHeight.SHORT_NAME, LongName = og.FigureHeight.LONG_NAME, Description = og.FigureHeight.DESCRIPTION, ValueName = "")]
        [Range(og.FigureHeight.MINIMUM, og.FigureHeight.MAXIMUM, ErrorMessage = og.FigureHeight.VALIDATION_ERROR_MESSAGE)]
        public int FigureHeight { get; set; }

        [Option(ShortName = og.XAxisMajorStep.SHORT_NAME, LongName = og.XAxisMajorStep.LONG_NAME, Description = og.XAxisMajorStep.DESCRIPTION, ValueName = "")]
        [Range(og.XAxisMajorStep.MINIMUM, og.XAxisMajorStep.MAXIMUM, ErrorMessage = og.XAxisMajorStep.VALIDATION_ERROR_MESSAGE)]
        public int XAxisMajorStep { get; set; }


        [Option(ShortName = io.DisplayAnnotationImpacts.SHORT_NAME, LongName = io.DisplayAnnotationImpacts.LONG_NAME,
            Description = io.DisplayAnnotationImpacts.DESCRIPTION, ValueName = "")]
        public string DisplayAnnotationImpacts { get; set; }

        [Option(ShortName = PolyploidQtlSeqCore.Share.ThreadNumber.SHORT_NAME, LongName = PolyploidQtlSeqCore.Share.ThreadNumber.LONG_NAME,
            Description = PolyploidQtlSeqCore.Share.ThreadNumber.DESCRIPTION, ValueName = "")]
        [Range(PolyploidQtlSeqCore.Share.ThreadNumber.MINIMUM, PolyploidQtlSeqCore.Share.ThreadNumber.MAXIMUM, ErrorMessage = PolyploidQtlSeqCore.Share.ThreadNumber.VALIDATION_ERROR_MESSAGE)]
        public int ThreadNumber { get; set; }

        [Option(ShortName = op.ParameterFileParser.SHORT_NAME, LongName = op.ParameterFileParser.LONG_NAME,
            Description = op.ParameterFileParser.DESCRIPTION, ValueName = "")]
        public string ParameterFile { get; set; }



        /// <summary>
        /// バージョン番号を取得する。
        /// </summary>
        /// <returns>Ver情報</returns>
        private static string GetVersion()
        {
            var version = typeof(PolyploidQtlSeqCommand).Assembly.GetCustomAttribute<AssemblyInformationalVersionAttribute>()?.InformationalVersion;

            return "Polyploid QTL-seq Ver " + version;
        }

        public override async Task<int> OnExecuteAsync(CommandLineApplication app)
        {
            int code;
            try
            {
                var pipeline = new QtlSeqPipeline(this, app.Options);
                code = await pipeline.RunAsync();
            }
            catch (Exception ex)
            {
                Console.Error.WriteLine(ex.Message);
                code = 1;
            }

            return code;
        }
    }
}
