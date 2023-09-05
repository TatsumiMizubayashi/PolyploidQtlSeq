using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeq.Options.Pipeline;
using PolyploidQtlSeq.Options.QtlAnalysis;
using PolyploidQtlSeqCore.Application.Pipeline;
using PolyploidQtlSeqCore.Options;
using System.Reflection;

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
    internal sealed class PolyploidQtlSeqCommand : CommandBase, IQtlSeqPipelineOptionValue
    {
        private static readonly string _parameterFileTitle = "QTL-Seq pipeline Parameter file";

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
            ChrSizeThreshold = ChrSizeThresholdOption.DEFAULT;
            AnalysisChrNames = "";
            MinMq = MinMappingQualityOption.DEFAULT;
            MinBq = MinBaseQualityOption.DEFAULT;
            AdjustMq = AdjustMappingQualityOption.DEFAULT;
            SnpEffMaxHeap = SnpEffMaxHeapSizeOption.DEFAULT;
            SnpEffConfigFile = "";
            SnpEffDatabaseName = "";
            OutputDir = "";
            Parent1MostAlleleRateThreshold = Parent1MostAlleleRateThresholdOption.DEFAULT;
            Parent2SnpIndexRange = Parent2SnpIndexRangeOption.DEFAULT;
            MinimumDepthThreshold = MinimumDepthThresholdOption.DEFAULT;
            MaxBulkSnpIndexThreshold = MaximumBulkSnpIndexThresholdOption.DEFAULT;
            Ploidy = PloidyOption.DEFAULT;
            Parent2PlexNumber = Parent2PlexNumberOption.DEFAULT;
            Bulk1Number = Bulk1NumberOption.DEFAULT;
            Bulk2Number = Bulk2NumberOption.DEFAULT;
            ReplicatesNumber = ReplicatesNumberOption.DEFAULT;
            WindowSize = WindowSizeOption.DEFAULT;
            StepSize = StepSizeOption.DEFAULT;
            FigureWidth = FigureWidthOption.DEFAULT;
            FigureHeight = FigureHeightOption.DEFAULT;
            XAxisMajorStep = XAxisMajorStepOption.DEFAULT;
            DisplayAnnotationImpacts = DisplayAnnotationImpactsOption.DEFAULT;
            ThreadNumber = ThreadNumberOption.DEFAULT;
            ParameterFilePath = "";
        }

        [Option(ShortName = ReferenceSequenceFileOption.SHORT_NAME, LongName = ReferenceSequenceFileOption.LONG_NAME,
            Description = ReferenceSequenceFileOption.DESCRIPTION, ValueName = "")]
        public string ReferenceSequence { get; set; }

        [Option(ShortName = Parent1DirectoryOption.SHORT_NAME, LongName = Parent1DirectoryOption.LONG_NAME,
            Description = Parent1DirectoryOption.DESCRIPTION, ValueName = "")]
        public string Parent1Dir { get; set; }

        [Option(ShortName = Parent2DirectoryOption.SHORT_NAME, LongName = Parent2DirectoryOption.LONG_NAME,
            Description = Parent2DirectoryOption.DESCRIPTION, ValueName = "")]
        public string Parent2Dir { get; set; }

        [Option(ShortName = Bulk1DirectoryOption.SHORT_NAME, LongName = Bulk1DirectoryOption.LONG_NAME,
            Description = Bulk1DirectoryOption.DESCRIPTION, ValueName = "")]
        public string Bulk1Dir { get; set; }

        [Option(ShortName = Bulk2DirectoryOption.SHORT_NAME, LongName = Bulk2DirectoryOption.LONG_NAME,
            Description = Bulk2DirectoryOption.DESCRIPTION, ValueName = "")]
        public string Bulk2Dir { get; set; }

        [Option(ShortName = OutputDirectoryOption.SHORT_NAME, LongName = OutputDirectoryOption.LONG_NAME,
            Description = OutputDirectoryOption.DESCRIPTION, ValueName = "")]
        public string OutputDir { get; set; }

        [Option(ShortName = ChrSizeThresholdOption.SHORT_NAME, LongName = ChrSizeThresholdOption.LONG_NAME,
            Description = ChrSizeThresholdOption.DESCRIPTION, ValueName = "")]
        public int ChrSizeThreshold { get; set; }

        [Option(ShortName = AnalysisChrNamesOption.SHORT_NAME, LongName = AnalysisChrNamesOption.LONG_NAME,
            Description = AnalysisChrNamesOption.DESCRIPTION, ValueName = "")]
        public string AnalysisChrNames { get; set; }

        [Option(ShortName = MinMappingQualityOption.SHORT_NAME, LongName = MinMappingQualityOption.LONG_NAME,
            Description = MinMappingQualityOption.DESCRIPTION, ValueName = "")]
        public int MinMq { get; set; }

        [Option(ShortName = MinBaseQualityOption.SHORT_NAME, LongName = MinBaseQualityOption.LONG_NAME,
            Description = MinBaseQualityOption.DESCRIPTION, ValueName = "")]
        public int MinBq { get; set; }

        [Option(ShortName = AdjustMappingQualityOption.SHORT_NAME, LongName = AdjustMappingQualityOption.LONG_NAME,
            Description = AdjustMappingQualityOption.DESCRIPTION, ValueName = "")]
        public int AdjustMq { get; set; }

        [Option(ShortName = SnpEffConfigFileOption.SHORT_NAME, LongName = SnpEffConfigFileOption.LONG_NAME,
            Description = SnpEffConfigFileOption.DESCRIPTION, ValueName = "")]
        public string SnpEffConfigFile { get; set; }

        [Option(ShortName = SnpEffDatabaseNameOption.SHORT_NAME, LongName = SnpEffDatabaseNameOption.LONG_NAME,
            Description = SnpEffDatabaseNameOption.DESCRIPTION, ValueName = "")]
        public string SnpEffDatabaseName { get; set; }

        [Option(ShortName = SnpEffMaxHeapSizeOption.SHORT_NAME, LongName = SnpEffMaxHeapSizeOption.LONG_NAME,
            Description = SnpEffMaxHeapSizeOption.DESCRIPTION, ValueName = "")]
        public int SnpEffMaxHeap { get; set; }


        [Option(ShortName = Parent1MostAlleleRateThresholdOption.SHORT_NAME, LongName = Parent1MostAlleleRateThresholdOption.LONG_NAME,
            Description = Parent1MostAlleleRateThresholdOption.DESCRIPTION, ValueName = "")]
        public double Parent1MostAlleleRateThreshold { get; set; }

        [Option(ShortName = Parent2SnpIndexRangeOption.SHORT_NAME, LongName = Parent2SnpIndexRangeOption.LONG_NAME,
            Description = Parent2SnpIndexRangeOption.DESCRIPTION, ValueName = "")]
        public string Parent2SnpIndexRange { get; set; }

        [Option(ShortName = MinimumDepthThresholdOption.SHORT_NAME, LongName = MinimumDepthThresholdOption.LONG_NAME,
            Description = MinimumDepthThresholdOption.DESCRIPTION, ValueName = "")]
        public int MinimumDepthThreshold { get; set; }

        [Option(ShortName = MaximumBulkSnpIndexThresholdOption.SHORT_NAME, LongName = MaximumBulkSnpIndexThresholdOption.LONG_NAME,
            Description = MaximumBulkSnpIndexThresholdOption.DESCRIPTION, ValueName = "")]
        public double MaxBulkSnpIndexThreshold { get; set; }


        [Option(ShortName = PloidyOption.SHORT_NAME, LongName = PloidyOption.LONG_NAME,
            Description = PloidyOption.DESCRIPTION, ValueName = "")]
        public int Ploidy { get; set; }

        [Option(ShortName = Parent2PlexNumberOption.SHORT_NAME, LongName = Parent2PlexNumberOption.LONG_NAME,
            Description = Parent2PlexNumberOption.DESCRIPTION, ValueName = "")]
        public int Parent2PlexNumber { get; set; }

        [Option(ShortName = Bulk1NumberOption.SHORT_NAME, LongName = Bulk1NumberOption.LONG_NAME,
            Description = Bulk1NumberOption.DESCRIPTION, ValueName = "")]
        public int Bulk1Number { get; set; }

        [Option(ShortName = Bulk2NumberOption.SHORT_NAME, LongName = Bulk2NumberOption.LONG_NAME,
            Description = Bulk2NumberOption.DESCRIPTION, ValueName = "")]
        public int Bulk2Number { get; set; }

        [Option(ShortName = ReplicatesNumberOption.SHORT_NAME, LongName = ReplicatesNumberOption.LONG_NAME,
            Description = ReplicatesNumberOption.DESCRIPTION, ValueName = "")]
        public int ReplicatesNumber { get; set; }

        [Option(ShortName = WindowSizeOption.SHORT_NAME, LongName = WindowSizeOption.LONG_NAME,
            Description = WindowSizeOption.DESCRIPTION, ValueName = "")]
        public int WindowSize { get; set; }

        [Option(ShortName = StepSizeOption.SHORT_NAME, LongName = StepSizeOption.LONG_NAME,
            Description = StepSizeOption.DESCRIPTION, ValueName = "")]
        public int StepSize { get; set; }

        [Option(ShortName = FigureWidthOption.SHORT_NAME, LongName = FigureWidthOption.LONG_NAME,
            Description = FigureWidthOption.DESCRIPTION, ValueName = "")]
        public int FigureWidth { get; set; }

        [Option(ShortName = FigureHeightOption.SHORT_NAME, LongName = FigureHeightOption.LONG_NAME, 
            Description = FigureHeightOption.DESCRIPTION, ValueName = "")]
        public int FigureHeight { get; set; }

        [Option(ShortName = XAxisMajorStepOption.SHORT_NAME, LongName = XAxisMajorStepOption.LONG_NAME, 
            Description = XAxisMajorStepOption.DESCRIPTION, ValueName = "")]
        public int XAxisMajorStep { get; set; }


        [Option(ShortName = DisplayAnnotationImpactsOption.SHORT_NAME, LongName = DisplayAnnotationImpactsOption.LONG_NAME,
            Description = DisplayAnnotationImpactsOption.DESCRIPTION, ValueName = "")]
        public string DisplayAnnotationImpacts { get; set; }

        [Option(ShortName = ThreadNumberOption.SHORT_NAME, LongName = ThreadNumberOption.LONG_NAME,
            Description = ThreadNumberOption.DESCRIPTION, ValueName = "")]
        public int ThreadNumber { get; set; }

        [Option(ShortName = ParameterFile.SHORT_NAME, LongName = ParameterFile.LONG_NAME,
            Description = ParameterFile.DESCRIPTION, ValueName = "")]
        public string ParameterFilePath { get; set; }

        public string InputVcf { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }


        public override async Task<int> OnExecuteAsync(CommandLineApplication app)
        {
            var options = QtlSeqPipelineOptions.Create(this);
            LoadParameterFile(ParameterFilePath, _parameterFileTitle, options, app);
            if (Validation(options)) Environment.Exit(1);

            int code;
            try
            {
                var pipeline = CreateQtlSeqPipeline(this);
                code = await pipeline.RunAsync();
            }
            catch (Exception ex)
            {
                Console.Error.WriteLine(ex.Message);
                code = 1;
            }
            finally
            {
                var logParamsFilePath = Path.Combine(OutputDir, "pipeline.params.txt");
                CreateParameterFile(logParamsFilePath, _parameterFileTitle, options);
            }

            return code;
        }

        private static QtlSeqPipeline CreateQtlSeqPipeline(IQtlSeqPipelineOptionValue optionValue)
        {
            var variantCallPipelineSettingValue = optionValue.CreateVariantCallPipelineSettingValue();
            var qtlSeqAnalysisSettingValue = optionValue.CreateQtlAnalysisScenarioSettingValue();

            return new QtlSeqPipeline(variantCallPipelineSettingValue, qtlSeqAnalysisSettingValue);
        }


        /// <summary>
        /// バージョン番号を取得する。
        /// </summary>
        /// <returns>Ver情報</returns>
        private static string GetVersion()
        {
            var version = typeof(PolyploidQtlSeqCommand).Assembly.GetCustomAttribute<AssemblyInformationalVersionAttribute>()?.InformationalVersion;

            return "Polyploid QTL-seq Ver " + version;
        }
    }
}
