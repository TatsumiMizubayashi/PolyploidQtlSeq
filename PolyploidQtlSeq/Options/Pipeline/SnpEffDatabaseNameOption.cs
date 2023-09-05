using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// SnpEff データベース名オプション
    /// </summary>
    internal class SnpEffDatabaseNameOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "sd";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "snpEffDatabase";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "SnpEff database name.";

        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// SnpEff データベース名オプション インスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public SnpEffDatabaseNameOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.SnpEffDatabaseName;

        protected override void SetValue(string value) => _optionValue.SnpEffDatabaseName = value;
    }
}
