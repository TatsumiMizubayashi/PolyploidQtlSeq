using PolyploidQtlSeqCore.Options;


namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// 出力ディレクトリオプション
    /// </summary>
    internal class OutputDirectoryOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "o";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "outputDir";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Output Directory.";

        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// 出力ディレクトリオプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public OutputDirectoryOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if (string.IsNullOrEmpty(_settingValue.OutputDir))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "Specify the output directory.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.OutputDir;

        protected override void SetValue(string value) => _settingValue.OutputDir = value;
    }
}
