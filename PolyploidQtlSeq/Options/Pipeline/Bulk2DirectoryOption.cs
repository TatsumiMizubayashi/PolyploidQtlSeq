using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// Bulk2ディレクトリオプション
    /// </summary>
    internal class Bulk2DirectoryOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "b2";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "bulk2";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Bulk2 Directory.";


        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// Bulk2ディレクトリオプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public Bulk2DirectoryOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (string.IsNullOrEmpty(_optionValue.Bulk2Dir))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "Specify a bulk2 Directory.");

            if (!Directory.Exists(_optionValue.Bulk2Dir))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"{_optionValue.Bulk2Dir} not found.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.Bulk2Dir;

        protected override void SetValue(string value) => _optionValue.Bulk2Dir = value;
    }
}
