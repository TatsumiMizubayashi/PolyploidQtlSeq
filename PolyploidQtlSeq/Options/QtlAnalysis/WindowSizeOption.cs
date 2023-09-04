using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// Window Sizeオプション
    /// </summary>
    internal class WindowSizeOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "w";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "window";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Window size (kbp) of the sliding window analysis.";

        /// <summary>
        /// window sizeの規定値
        /// </summary>
        public const int DEFAULT = 100;

        /// <summary>
        /// window sizeの最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// window sizeの最大値
        /// </summary>
        private const int MAXIMUM = 100000;

        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// Window Sizeオプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public WindowSizeOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _settingValue.WindowSize && _settingValue.WindowSize <= MAXIMUM)
                return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"Window size (kbp) should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.WindowSize.ToString();

        protected override void SetValue(string value) => _settingValue.WindowSize = int.Parse(value);
    }
}
