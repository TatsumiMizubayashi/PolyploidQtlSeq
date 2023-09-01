using PolyploidQtlSeqCore.Application.QualityControl;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QualityControl
{
    /// <summary>
    /// 入力RawFastqディレクトリオプション
    /// </summary>
    internal class InputRawFastqDirectoryOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "i";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "inputDir";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Raw fastq directory.";

        private readonly IFastpQualityControlSettingValue _optionValue;

        /// <summary>
        /// 入力RawFastqディレクトリオプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">Fastq QCオプション値</param>
        public InputRawFastqDirectoryOption(IFastpQualityControlSettingValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (string.IsNullOrEmpty(_optionValue.InputDir))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "Specify the input raw Fastq directory.");

            if (!Directory.Exists(_optionValue.InputDir))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"{_optionValue.InputDir} not found.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.InputDir;

        protected override void SetValue(string value) => _optionValue.InputDir = value;
    }
}
