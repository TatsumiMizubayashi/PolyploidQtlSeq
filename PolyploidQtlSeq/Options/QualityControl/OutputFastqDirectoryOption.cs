using PolyploidQtlSeqCore.Application.QualityControl;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QualityControl
{
    /// <summary>
    /// 出力Fastqディレクトリオプション
    /// </summary>
    internal class OutputFastqDirectoryOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "o";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "outputDir";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Output directory.";

        private readonly IFastpQualityControlOptionValue _optionValue;

        /// <summary>
        /// 出力Fastqディレクトリオプション インスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">Fastq QCオプション値</param>
        public OutputFastqDirectoryOption(IFastpQualityControlOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (string.IsNullOrEmpty(_optionValue.OutputDir))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "Specify the output directory.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.OutputDir;

        protected override void SetValue(string value) => _optionValue.OutputDir = value;
    }
}
