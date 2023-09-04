using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// StepSizeオプション
    /// </summary>
    internal class StepSizeOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "s";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "step";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Step size (kbp) of the sliding window analysis.";

        /// <summary>
        /// setp sizeの規定値
        /// </summary>
        public const int DEFAULT = 20;

        /// <summary>
        /// step sizeの最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// step sizeの最大値
        /// </summary>
        private const int MAXIMUM = 100000;

        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// StepSizeオプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public StepSizeOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if(_settingValue.StepSize < MINIMUM || _settingValue.StepSize > MAXIMUM)
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"Step size (kbp) should be an integer between {MINIMUM} and {MAXIMUM}.");

            if(_settingValue.StepSize > _settingValue.WindowSize)
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"The step size should be less than or equal to the window size.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.StepSize.ToString();

        protected override void SetValue(string value) => _settingValue.StepSize = int.Parse(value);
    }
}
