using McMaster.Extensions.CommandLineUtils;
using PolyploidQtlSeq.Options.QtlAnalysis;
using PolyploidQtlSeqCore.Application.QtlAnalysis;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq
{
    /// <summary>
    /// QTL-Seq解析を行うコマンド
    /// </summary>
    [Command("qtl", Description = "QTL-Seq analysis.")]
    internal sealed class QtlAnalysisCommand : CommandBase, IQtlSeqAnalysisOptionValue
    {
        private static readonly string _parameterFileTitle = "QTL-Seq analysis Parameter file";

        /// <summary>
        /// QTL解析コマンドを作成する。
        /// </summary>
        public QtlAnalysisCommand()
        {
            InputVcf = "";
            OutputDir = "";
            Parent1MostAlleleRateThreshold = Parent1MostAlleleRateThresholdOption.DEFAULT;
            Parent2SnpIndexRange = Parent2SnpIndexRangeOption.GetDefault();
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

        [Option(ShortName = InputVcfFileOption.SHORT_NAME, LongName = InputVcfFileOption.LONG_NAME,
            Description = InputVcfFileOption.DESCRIPTION, ValueName = "")]
        public string InputVcf { get; set; }

        [Option(ShortName = OutputDirectoryOption.SHORT_NAME, LongName = OutputDirectoryOption.LONG_NAME,
            Description = OutputDirectoryOption.DESCRIPTION, ValueName = "")]
        public string OutputDir { get; set; }

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


        public override async Task<int> OnExecuteAsync(CommandLineApplication app)
        {
            var options = QtlSeqAnalysisOptions.Create(this);
            LoadParameterFile(ParameterFilePath, _parameterFileTitle, options, app);
            if (Validation(options)) Environment.Exit(1);

            int code;
            try
            {
                code = await RunAsync();
            }
            catch(Exception ex)
            {
                Console.Error.WriteLine(ex.Message);
                code = 1;
            }
            finally
            {
                var logParamsFilePath = Path.Combine(OutputDir, "qtl.params.txt");
                CreateParameterFile(logParamsFilePath, _parameterFileTitle, options);
            }

            return code;
        }

        // OnExecuteが非同期なのでTaskを使用して処理する
        private async ValueTask<int> RunAsync()
        {
            return await Task.Run(() =>
            {
                var analyzer = CreateQtlSeqAnalysis(this);
                return analyzer.Run(InputVcf);
            });
        }

        /// <summary>
        /// QtlSeq解析インスタンスを作成する。
        /// </summary>
        /// <param name="optioValue">設定値</param>
        /// <returns>QtlSeq解析インスタンス</returns>
        private static QtlSeqAnalysis CreateQtlSeqAnalysis(IQtlSeqAnalysisOptionValue optioValue)
        {
            var setting = optioValue.CreateQtlAnalysisScenarioSettingValue();
            return new QtlSeqAnalysis(setting);
        }

    }
}
