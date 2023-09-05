using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// Bulk1ディレクトリオプション
    /// </summary>
    internal class Bulk1DirectoryOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "b1";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "bulk1";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Bulk1 Directory.";

        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// Bulk1ディレクトリオプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public Bulk1DirectoryOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (string.IsNullOrEmpty(_optionValue.Bulk1Dir))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "Specify a bulk1 Directory.");

            if (!Directory.Exists(_optionValue.Bulk1Dir))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"{_optionValue.Bulk1Dir} not found.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.Bulk1Dir;

        protected override void SetValue(string value) => _optionValue.Bulk1Dir = value;
    }
}
