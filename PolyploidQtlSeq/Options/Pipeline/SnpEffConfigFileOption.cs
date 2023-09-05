using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.Pipeline
{
    /// <summary>
    /// SnpEff.configファイルオプション
    /// </summary>
    internal class SnpEffConfigFileOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "sc";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "snpEffConfig";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "snpEff.config file. Not required if snpEff default config file is used.";

        private readonly IQtlSeqPipelineOptionValue _optionValue;

        /// <summary>
        /// SnpEff.configファイルオプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">QtlSeqパイプラインオプション値</param>
        public SnpEffConfigFileOption(IQtlSeqPipelineOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (string.IsNullOrEmpty(_optionValue.SnpEffConfigFile)) return new DataValidationResult();
            if (File.Exists(_optionValue.SnpEffConfigFile)) return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME, $"{_optionValue.SnpEffConfigFile} not found."); ;
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.SnpEffConfigFile;

        protected override void SetValue(string value) => _optionValue.SnpEffConfigFile = value;
    }
}
