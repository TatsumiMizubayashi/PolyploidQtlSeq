using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// 染色体長しきい値オプション
    /// </summary>
    internal class ChrSizeThresholdOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "cs";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "chrSize";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Threshold for length of chromosomes to be analyzed. Chromosomes with a length more than this value are analyzed.";

        /// <summary>
        /// しきい値の規定値
        /// </summary>
        public const int DEFAULT = 10_000_000;

        /// <summary>
        /// しきい値の最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// しきい値の最大値
        /// </summary>
        private const int MAXIMUM = 100_000_000;


        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// 染色体長しきい値オプション インスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public ChrSizeThresholdOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (_optionValue.ChrSizeThreshold < MINIMUM || _optionValue.ChrSizeThreshold > MAXIMUM)
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"Threshold for length of chromosomes to be analyzed should be an integer between {MINIMUM} and {MAXIMUM}.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.ChrSizeThreshold.ToString();

        protected override void SetValue(string value) => _optionValue.ChrSizeThreshold = int.Parse(value);
    }
}
