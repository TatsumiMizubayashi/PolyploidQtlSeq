using PolyploidQtlSeqCore.Application.QtlAnalysis;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// Bulk2個体数オプション
    /// </summary>
    internal class Bulk2NumberOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "n2";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "NBulk2";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Number of Individuals in bulk2.";

        /// <summary>
        /// 個体数の規定値
        /// </summary>
        public const int DEFAULT = 20;

        /// <summary>
        /// 個体数の最小値
        /// </summary>
        private const int MINIMUM = 2;

        /// <summary>
        /// 個体数の最大値
        /// </summary>
        private const int MAXIMUM = 1000;

        private readonly IQtlSeqAnalysisSettingValue _settingValue;

        /// <summary>
        /// Bulk2個体数オプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public Bulk2NumberOption(IQtlSeqAnalysisSettingValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _settingValue.Bulk2Number && _settingValue.Bulk2Number <= MAXIMUM)
                return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"Number of Individuals in bulk2 should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.Bulk2Number.ToString();

        protected override void SetValue(string value) => _settingValue.Bulk2Number = int.Parse(value);
    }
}
