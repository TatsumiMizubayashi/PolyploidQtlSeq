using PolyploidQtlSeqCore.Application.QtlAnalysis;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// 入力VCFファイルオプション
    /// </summary>
    internal class InputVcfFileOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "i";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "inputVcf";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Input VCF file.";

        private readonly IQtlSeqAnalysisSettingValue _settingValue;

        /// <summary>
        /// 入力VCFファイルオプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public InputVcfFileOption(IQtlSeqAnalysisSettingValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if (string.IsNullOrEmpty(_settingValue.InputVcf))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "Specify a VCF file.");

            if (!File.Exists(_settingValue.InputVcf))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"{_settingValue.InputVcf} not found.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.InputVcf;

        protected override void SetValue(string value) => _settingValue.InputVcf = value;
    }
}
