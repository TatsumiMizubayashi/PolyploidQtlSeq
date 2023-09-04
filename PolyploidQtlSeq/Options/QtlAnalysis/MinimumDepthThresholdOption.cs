using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// 最小Depthしきい値オプション
    /// </summary>
    internal class MinimumDepthThresholdOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "md";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "minDepth";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Minimum Depth threshold. The variants with even one sample below this threshold are excluded for QTL analysis.";

        /// <summary>
        /// 最低Depthしきい値の規定値
        /// </summary>
        public const int DEFAULT = 40;

        /// <summary>
        /// 最低Depthしきい値の最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 最低Depthしきい値の最大値
        /// </summary>
        private const int MAXIMUM = 10000;

        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// 最小Depthしきい値オプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public MinimumDepthThresholdOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _settingValue.MinimumDepthThreshold && _settingValue.MinimumDepthThreshold <= MAXIMUM)
                return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME, 
                $"Minimum Depth threshold should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.MinimumDepthThreshold.ToString();

        protected override void SetValue(string value) => _settingValue.MinimumDepthThreshold = int.Parse(value);
    }
}
