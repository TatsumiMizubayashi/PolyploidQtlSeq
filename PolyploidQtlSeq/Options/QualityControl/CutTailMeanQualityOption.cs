using PolyploidQtlSeqCore.Application.QualityControl;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QualityControl
{
    /// <summary>
    /// 3'末端トリム 平均クオリティオプション
    /// </summary>
    internal class CutTailMeanQualityOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "Q";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "cutTailMeanQuality";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Threshold of average quality for trimming at 3' end. cut_tail_mean_quality in fastp.";

        /// <summary>
        /// 3'末端トリム時の平均クオリティの規定値
        /// </summary>
        public const int DEFAULT = 20;

        /// <summary>
        /// 3'末端トリム時の平均クオリティの最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 3'末端トリム時の平均クオリティの最大値
        /// </summary>
        private const int MAXIMUM = 30;

        private readonly IFastpQualityControlSettingValue _optionValue;

        /// <summary>
        /// 3'末端トリム 平均クオリティオプション インスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">Fastq QCオプション値</param>
        public CutTailMeanQualityOption(IFastpQualityControlSettingValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _optionValue.CutTailMeanQuality && _optionValue.CutTailMeanQuality <= MAXIMUM) return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"Threshold of average quality for trimming at 3' end should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.CutTailMeanQuality.ToString();

        protected override void SetValue(string value) => _optionValue.CutTailMeanQuality = int.Parse(value);
    }
}
