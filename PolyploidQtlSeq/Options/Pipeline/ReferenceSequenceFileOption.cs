using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// リファレンスシークエンスファイルオプション
    /// </summary>
    internal class ReferenceSequenceFileOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "r";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "refSeq";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Reference sequence file.";

        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// リファレンスシークエンスファイルオプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public ReferenceSequenceFileOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (string.IsNullOrEmpty(_optionValue.ReferenceSequence))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "Specify a reference sequence file.");

            if (!File.Exists(_optionValue.ReferenceSequence))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"{_optionValue.ReferenceSequence} not found.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.ReferenceSequence;

        protected override void SetValue(string value) => _optionValue.ReferenceSequence = value;
    }
}
