using PolyploidQtlSeqCore.Application.QualityControl;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QualityControl
{
    /// <summary>
    /// Fastp QCオプション
    /// </summary>
    internal class FastpQualityControlOptions : OptionCollection
    {
        /// <summary>
        /// Fastp QCオプションインスタンスを作成する。
        /// </summary>
        /// <param name="options">オプション</param>
        private FastpQualityControlOptions(Option[] options) : base(options)
        {
        }

        /// <summary>
        /// Fastp QCオプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">Fastp QCオプション値</param>
        /// <returns>Fastp QCオプション</returns>
        public static FastpQualityControlOptions Create(IFastpQualityControlSettingValue optionValue)
        {
            var options = new Option[]
            {
                new InputRawFastqDirectoryOption(optionValue),
                new OutputFastqDirectoryOption(optionValue),
                new ReadLengthRequiredOption(optionValue),
                new NBaseLimitOption(optionValue),
                new BaseQualityOption(optionValue),
                new CutTailMeanQualityOption(optionValue),
                new CutTailWindowSizeOption(optionValue),
                new ThreadNumberOption(optionValue),
            };

            return new FastpQualityControlOptions(options);
        }
    }
}
