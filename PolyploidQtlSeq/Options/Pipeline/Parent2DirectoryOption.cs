using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// P2 ディレクトリオプション
    /// </summary>
    internal class Parent2DirectoryOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "p2";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "parent2";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Parent2 Directory.";

        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// P2ディレクトリオプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public Parent2DirectoryOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (string.IsNullOrEmpty(_optionValue.Parent2Dir))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "Specify a parent2 Directory.");

            if (!Directory.Exists(_optionValue.Parent2Dir))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"{_optionValue.Parent2Dir} not found.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.Parent2Dir;

        protected override void SetValue(string value) => _optionValue.Parent2Dir = value;
    }
}
