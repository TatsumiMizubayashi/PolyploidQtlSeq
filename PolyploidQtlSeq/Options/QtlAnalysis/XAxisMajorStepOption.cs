using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// X軸目盛り間隔オプション
    /// </summary>
    internal class XAxisMajorStepOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "xs";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "xStep";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "X-axis scale interval (Mbp) of the graphs.";

        /// <summary>
        /// 目盛り間隔の規定値
        /// </summary>
        public const int DEFAULT = 5;

        /// <summary>
        /// 目盛り間隔の最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 目盛り間隔の最大値
        /// </summary>
        private const int MAXIMUM = 50;

        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// X軸目盛り間隔オプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public XAxisMajorStepOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _settingValue.XAxisMajorStep && _settingValue.XAxisMajorStep <= MAXIMUM)
                return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"X-axis scale interval (Mbp) of the graphs should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.XAxisMajorStep.ToString();

        protected override void SetValue(string value) => _settingValue.XAxisMajorStep = int.Parse(value);
    }
}
