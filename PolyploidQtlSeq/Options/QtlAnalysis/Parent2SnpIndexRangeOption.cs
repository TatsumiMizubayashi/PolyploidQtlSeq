using PolyploidQtlSeqCore.Application.QtlAnalysis;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// P2 SNP-index範囲オプション
    /// </summary>
    internal class Parent2SnpIndexRangeOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "p2r";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "p2SnpIndexRange";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "SNP-index range for parent2.";

        /// <summary>
        /// Parent2 SNP-index範囲の規定値
        /// </summary>
        public const string DEFAULT = "0.15-0.375";

        /// <summary>
        /// SNP-indexの最小値
        /// </summary>
        private const double MINIMUM = 0;

        /// <summary>
        /// SNP-indexの最大値
        /// </summary>
        private const double MAXIMUM = 1.0;

        private static readonly char _delimiter = '-';

        private readonly IQtlSeqAnalysisSettingValue _settingValue;

        /// <summary>
        /// P2 SNP-index範囲オプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public Parent2SnpIndexRangeOption(IQtlSeqAnalysisSettingValue settingValue)
        {
            _settingValue = settingValue;
        }   

        public override DataValidationResult Validation()
        {
            if(string.IsNullOrEmpty(_settingValue.Parent2SnpIndexRange)) 
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "Specify by lower limit - upper limit.");

            var values = _settingValue.Parent2SnpIndexRange.Split(_delimiter);
            if(values.Length != 2) return new DataValidationResult(SHORT_NAME, LONG_NAME, "Specify by lower limit - upper limit.");

            if (!double.TryParse(values[0], out var lower))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "The lower limit connot bo converted to a numerical value.");
            if (!double.TryParse(values[1], out var upper))
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "The upper limit cannot be converted to a numerical value.");

            if (lower < MINIMUM || lower > MAXIMUM) 
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "The lower limit should be specified in the range of 0.0 to 1.0.");
            if (upper < MINIMUM || upper > MAXIMUM)
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "The upper limit should be specified in the range of 0.0 to 1.0.");

            if (lower >= upper) return new DataValidationResult(SHORT_NAME, LONG_NAME, "The upper limit should be greater than the lower limit.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.Parent2SnpIndexRange;

        protected override void SetValue(string value) => _settingValue.Parent2SnpIndexRange = value;

    }
}
