using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// Adjust MQオプション
    /// </summary>
    internal class AdjustMappingQualityOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "C";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "adjustMQ";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Value for adjust mapping quality at variant detection in bcftools mpileup. \r\nSpecify 0, to disable this function.\r\n";

        /// <summary>
        /// adjust MQの規定値
        /// </summary>
        public const int DEFAULT = 60;

        /// <summary>
        /// adjust MQの最小値
        /// </summary>
        private const int MINIMUM = 0;

        /// <summary>
        /// adjust MQの最大値
        /// </summary>
        private const int MAXIMUM = 1000;

        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// Adjust MQオプション インスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public AdjustMappingQualityOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (_optionValue.AdjustMq < MINIMUM || _optionValue.AdjustMq > MAXIMUM)
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"AdjustMQ should be an integer between {MINIMUM} and {MAXIMUM}.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.AdjustMq.ToString();

        protected override void SetValue(string value) => _optionValue.AdjustMq = int.Parse(value);
    }
}
