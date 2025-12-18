using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// P1最多アレル割合しきい値オプション
    /// </summary>
    internal class Parent1MostAlleleRateThresholdOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "p1r";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "p1MostAlleleRate";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Most allele frequency for Parent1. Variants exceeding this threshold is considered homozygous.";

        /// <summary>
        /// Parent1 MaxAllelRateしきい値の規定値
        /// </summary>
        public const double DEFAULT = 0.99;

        /// <summary>
        /// Parent1 MaxAllelRateしきい値の最小値
        /// </summary>
        private const double MINIMUM = 0;

        /// <summary>
        /// Parent1 MaxAllelRateしきい値の最大値
        /// </summary>
        private const double MAXIMUM = 1.0;

        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// P1最多アレル割合しきい値オプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public Parent1MostAlleleRateThresholdOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _settingValue.Parent1MostAlleleRateThreshold && _settingValue.Parent1MostAlleleRateThreshold <= MAXIMUM)
                return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"Most allele frequency for Parent1 should be an integer between {MINIMUM:F1} and {MAXIMUM:F1}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.Parent1MostAlleleRateThreshold.ToString();

        protected override void SetValue(string value) => _settingValue.Parent1MostAlleleRateThreshold = double.Parse(value);
    }
}
