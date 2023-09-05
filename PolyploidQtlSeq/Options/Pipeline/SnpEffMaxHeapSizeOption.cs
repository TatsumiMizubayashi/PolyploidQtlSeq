using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// SnpEff最大ヒープサイズオプション
    /// </summary>
    internal class SnpEffMaxHeapSizeOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "sm";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "snpEffMaxHeap";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "SnpEff maximum heap size (GB).";

        /// <summary>
        /// max heapの規定値
        /// </summary>
        public const int DEFAULT = 6;

        /// <summary>
        /// max heapの最小値
        /// </summary>
        private const int MINIMUM = 2;

        /// <summary>
        /// max heapの最大値
        /// </summary>
        private const int MAXIMUM = 100;

        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// SnpEff最大ヒープサイズオプション インスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public SnpEffMaxHeapSizeOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (_optionValue.SnpEffMaxHeap < MINIMUM || _optionValue.SnpEffMaxHeap > MAXIMUM)
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"SnpEff maximum heap size (GB) should be an integer between {MINIMUM} and {MAXIMUM}.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.SnpEffMaxHeap.ToString();

        protected override void SetValue(string value) => _optionValue.SnpEffMaxHeap = int.Parse(value);
    }
}
