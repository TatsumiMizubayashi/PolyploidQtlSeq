using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// P1ディレクトリオプション
    /// </summary>
    internal class Parent1DirectoryOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "p1";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "parent1";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Parent1 Directory.";

        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// P1ディレクトリオプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public Parent1DirectoryOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (string.IsNullOrEmpty(_optionValue.Parent1Dir))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "Specify a parent1 Directory.");

            if (!Directory.Exists(_optionValue.Parent1Dir))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"{_optionValue.Parent1Dir} not found.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.Parent1Dir;

        protected override void SetValue(string value) => _optionValue.Parent1Dir = value;
    }
}
