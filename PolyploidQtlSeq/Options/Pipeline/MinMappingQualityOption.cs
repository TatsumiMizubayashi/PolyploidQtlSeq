using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// 最小Mappingクオリティオプション
    /// </summary>
    internal class MinMappingQualityOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "q";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "minMQ";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Minimum mapping quality at variant detection in bcftools mpileup.";

        /// <summary>
        /// min-MQの規定値
        /// </summary>
        public const int DEFAULT = 40;

        /// <summary>
        /// min-MQの最小値
        /// </summary>
        private const int MINIMUM = 0;

        /// <summary>
        /// min-MQの最大値
        /// </summary>
        private const int MAXIMUM = 60;

        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// 最小Mappingクオリティオプション インスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public MinMappingQualityOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (_optionValue.MinMq < MINIMUM || _optionValue.MinMq > MAXIMUM)
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"Minimum mapping quality should be an integer between {MINIMUM} and {MAXIMUM}.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.MinMq.ToString();

        protected override void SetValue(string value) => _optionValue.MinMq = int.Parse(value);
    }
}
