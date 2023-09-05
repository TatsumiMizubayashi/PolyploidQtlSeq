using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// 最小塩基クオリティオプション
    /// </summary>
    internal class MinBaseQualityOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "Q";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "minBQ";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Minimum base quality at variant detection in bcftools mpileup.";

        /// <summary>
        /// min-BQの規定値
        /// </summary>
        public const int DEFAULT = 13;

        /// <summary>
        /// min-BQの最小値
        /// </summary>
        private const int MINIMUM = 0;

        /// <summary>
        /// min-BQの最大値
        /// </summary>
        private const int MAXIMUM = 60;

        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// 最小塩基クオリティオプション インスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public MinBaseQualityOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (_optionValue.MinBq < MINIMUM || _optionValue.MinBq > MAXIMUM)
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"Minimum base quality should be an integer between {MINIMUM} and {MAXIMUM}.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.MinBq.ToString();

        protected override void SetValue(string value) => _optionValue.MinBq = int.Parse(value);
    }
}
