using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// 解析染色体名オプション
    /// </summary>
    internal class AnalysisChrNamesOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "cn";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "chrNames";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Specify the chromosome name to be analyzed. If there are more than one, separate them with commas. ";

        private static readonly char[] _delimiter = [','];

        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// 解析染色体名オプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public AnalysisChrNamesOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (string.IsNullOrEmpty(_optionValue.AnalysisChrNames)) return new DataValidationResult();

            var chrNames = _optionValue.AnalysisChrNames.Split(_delimiter);
            var uniqChrNames = chrNames.Distinct().ToArray();
            if (chrNames.Length == uniqChrNames.Length) return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME, "There is a duplicate chromosome name.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.AnalysisChrNames;

        protected override void SetValue(string value) => _optionValue.AnalysisChrNames = value;
    }
}
