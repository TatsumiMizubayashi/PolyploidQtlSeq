using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// 表示アノテーションインパクトオプション
    /// </summary>
    internal class DisplayAnnotationImpactsOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "di";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "displayImpacts";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Annotation Impact to be included in the SNP-index file. Separate multiple items with commas.";

        /// <summary>
        /// DisplayImpactsの規定値
        /// </summary>
        public const string DEFAULT = "HIGH,MODERATE";

        private static readonly char _delimiter = ',';

        private static readonly IReadOnlyDictionary<string, string> _impactDictionary = new Dictionary<string, string>()
        {
            ["HIGH"] = "",
            ["MODERATE"] = "",
            ["LOW"] = "",
            ["MODIFIER"] = ""
        };

        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// 表示アノテーションインパクトオプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public DisplayAnnotationImpactsOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if (string.IsNullOrEmpty(_settingValue.DisplayAnnotationImpacts)) return new DataValidationResult();

            var errorImpactValues = _settingValue.DisplayAnnotationImpacts.Split(_delimiter)
                .Select(x => x.ToUpper())
                .Where(x => !_impactDictionary.ContainsKey(x))
                .ToArray();
            if (errorImpactValues.Length == 0) return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME, "Select annotation impact from HIGH, MODERATE, LOW, or MODIFIER.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.DisplayAnnotationImpacts;

        protected override void SetValue(string value) => _settingValue.DisplayAnnotationImpacts = value;
    }
}
