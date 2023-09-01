using PolyploidQtlSeqCore.Application.QualityControl;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QualityControl
{
    /// <summary>
    /// 許容するN塩基数オプション
    /// </summary>
    internal class NBaseLimitOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "n";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "nBaseLimit";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Maximum number of N bases. n_base_limit in fastp.";

        /// <summary>
        /// N塩基数の規定値
        /// </summary>
        public const int DEFAULT = 5;

        /// <summary>
        /// N塩基数の最小値
        /// </summary>
        private const int MINIMUM = 0;

        /// <summary>
        /// N塩基数の最大値
        /// </summary>
        private const int MAXIMUM = 30;

        private readonly IFastpQualityControlSettingValue _optionValue;

        /// <summary>
        /// 許容するN塩基数オプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">Fastq QCオプション値</param>
        public NBaseLimitOption(IFastpQualityControlSettingValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _optionValue.NBaseLimit && _optionValue.NBaseLimit <= MAXIMUM) return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"The upper limit of the number of N bases should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.NBaseLimit.ToString();

        protected override void SetValue(string value) => _optionValue.NBaseLimit = int.Parse(value);
    }
}
