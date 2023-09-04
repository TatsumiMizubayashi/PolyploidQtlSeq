using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// 図表高オプション
    /// </summary>
    internal class FigureHeightOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "fh";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "figHeight";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Height (pixel) of the graph images.";

        /// <summary>
        /// グラフ高さの規定値
        /// </summary>
        public const int DEFAULT = 300;

        /// <summary>
        /// グラフ高さの最小値
        /// </summary>
        private const int MINIMUM = 100;

        /// <summary>
        /// グラフ高さの最大値
        /// </summary>
        private const int MAXIMUM = 2000;


        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// 図表高オプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public FigureHeightOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }   

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _settingValue.FigureHeight && _settingValue.FigureHeight <= MAXIMUM)
                return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"Height (pixel) of the graph images should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.FigureHeight.ToString();

        protected override void SetValue(string value) => _settingValue.FigureHeight = int.Parse(value);
    }
}
