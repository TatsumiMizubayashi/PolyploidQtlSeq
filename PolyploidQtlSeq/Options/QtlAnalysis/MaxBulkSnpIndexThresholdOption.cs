using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// 最大Bulk SNP-indexしきい値オプション
    /// </summary>
    internal class MaxBulkSnpIndexThresholdOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "mb";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "maxBulkSnpIndex";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Maximum threshold for SNP-index for the Bulk samples. Variants with a SNP-index exceeding this value are excluded.";

        /// <summary>
        /// 最大Bulk SNP-indexのしきい値の規定値
        /// </summary>
        public const double DEFAULT = 1.0;

        /// <summary>
        /// 最大Bulk SNP-indexのしきい値の最小値
        /// </summary>
        private const double MINIMUM = 0;

        /// <summary>
        /// 最大Bulk SNP-indexのしきい値の最大値
        /// </summary>
        private const double MAXIMUM = 1.0;

        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// 最大Bulk SNP-indexしきい値オプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public MaxBulkSnpIndexThresholdOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _settingValue.MaxBulkSnpIndexThreshold && _settingValue.MaxBulkSnpIndexThreshold <= MAXIMUM)
                return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"Maximum threshold for SNP-index for the Bulk samples should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.MaxBulkSnpIndexThreshold.ToString();

        protected override void SetValue(string value) => _settingValue.MaxBulkSnpIndexThreshold = double.Parse(value);
    }
}
