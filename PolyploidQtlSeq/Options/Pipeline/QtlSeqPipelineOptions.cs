using PolyploidQtlSeq.Options.QtlAnalysis;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// QTL-Seqパイプラインオプション
    /// </summary>
    internal class QtlSeqPipelineOptions : OptionCollection
    {
        /// <summary>
        /// QTL-Seqパイプラインオプションインスタンスを作成する。
        /// </summary>
        /// <param name="options">オプション</param>
        private QtlSeqPipelineOptions(Option[] options) : base(options)
        {
        }

        /// <summary>
        /// QTL-Seqパイプラインオプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QTL-Seqパイプラインオプション値</param>
        /// <returns>QTL-Seqパイプラインオプション</returns>
        public static QtlSeqPipelineOptions Create(IQtlSeqPipelineOptionValue optionValue)
        {
            var options = new Option[]
            {
                new ReferenceSequenceFileOption(optionValue),
                new Parent1DirectoryOption(optionValue),
                new Parent2DirectoryOption(optionValue),
                new Bulk1DirectoryOption(optionValue),
                new Bulk2DirectoryOption(optionValue),
                new OutputDirectoryOption(optionValue),
                new ChrSizeThresholdOption(optionValue),
                new AnalysisChrNamesOption(optionValue),
                new MinMappingQualityOption(optionValue),
                new MinBaseQualityOption(optionValue),
                new AdjustMappingQualityOption(optionValue),
                new SnpEffConfigFileOption(optionValue),
                new SnpEffDatabaseNameOption(optionValue),
                new SnpEffMaxHeapSizeOption(optionValue),
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

            return new QtlSeqPipelineOptions(options);
        }
    }
}
