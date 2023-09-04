using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// 図表幅オプション
    /// </summary>
    internal class FigureWidthOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "fw";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "figWidth";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Width (pixel) of the graph images.";

        /// <summary>
        /// グラフ幅の規定値
        /// </summary>
        public const int DEFAULT = 1200;

        /// <summary>
        /// グラフ幅の最小値
        /// </summary>
        private const int MINIMUM = 300;

        /// <summary>
        /// グラフ幅の最大値
        /// </summary>
        private const int MAXIMUM = 5000;

        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// 図表幅オプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public FigureWidthOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _settingValue.FigureWidth && _settingValue.FigureWidth <= MAXIMUM)
                return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"Width (pixel) of the graph images should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.FigureWidth.ToString();

        protected override void SetValue(string value) => _settingValue.FigureWidth = int.Parse(value);
    }
}
