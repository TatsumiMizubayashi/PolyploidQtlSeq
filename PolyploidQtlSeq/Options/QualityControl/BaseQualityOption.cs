using PolyploidQtlSeqCore.Application.QualityControl;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QualityControl
{
    /// <summary>
    /// 塩基クオリティオプション
    /// </summary>
    internal class BaseQualityOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "q";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "quality";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Threshold for base quality. qualified_quality_phred in fastp.";

        /// <summary>
        /// 塩基クオリティの規定値
        /// </summary>
        public const int DEFAULT = 15;

        /// <summary>
        /// 塩基クオリティの最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 塩基クオリティの最大値
        /// </summary>
        private const int MAXIMUM = 50;

        private readonly IFastpQualityControlOptionValue _optionValue;

        /// <summary>
        /// 塩基クオリティオプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">Fastq QCオプション値</param>
        public BaseQualityOption(IFastpQualityControlOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _optionValue.BaseQuality && _optionValue.BaseQuality <= MAXIMUM) return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"Base quality should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.BaseQuality.ToString();

        protected override void SetValue(string value) => _optionValue.BaseQuality = int.Parse(value);
    }
}
