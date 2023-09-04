using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// QTL-seq解析オプション
    /// </summary>
    internal class QtlSeqAnalysisOptions : OptionCollection
    {
        /// <summary>
        /// QTL-seq解析オプションインスタンスを作成する。
        /// </summary>
        /// <param name="options">オプション</param>
        private QtlSeqAnalysisOptions(Option[] options) : base(options)
        {
        }

        /// <summary>
        /// QTL-seq解析オプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QTL-Seq解析オプション値</param>
        /// <returns>QTL-seq解析オプション</returns>
        public static QtlSeqAnalysisOptions Create(IQtlSeqAnalysisOptionValue optionValue)
        {
            var options = new Option[]
            {
                new InputVcfFileOption(optionValue),
                new OutputDirectoryOption(optionValue),
                new Parent1MostAlleleRateThresholdOption(optionValue),
                new Parent2SnpIndexRangeOption(optionValue),
                new MinimumDepthThresholdOption(optionValue),
                new MaximumBulkSnpIndexThresholdOption(optionValue),
                new PloidyOption(optionValue),
                new Parent2PlexNumberOption(optionValue),
                new Bulk1NumberOption(optionValue),
                new Bulk2NumberOption(optionValue),
                new ReplicatesNumberOption(optionValue),
                new WindowSizeOption(optionValue),
                new StepSizeOption(optionValue),
                new FigureWidthOption(optionValue),
                new FigureHeightOption(optionValue),
                new XAxisMajorStepOption(optionValue),
                new DisplayAnnotationImpactsOption(optionValue),
                new ThreadNumberOption(optionValue),
            };

            return new QtlSeqAnalysisOptions(options);
        }
    }
}
