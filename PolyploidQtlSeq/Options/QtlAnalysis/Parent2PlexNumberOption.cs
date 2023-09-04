using PolyploidQtlSeqCore.Application.QtlAnalysis;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// P2 Plex数オプション
    /// </summary>
    internal class Parent2PlexNumberOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "np";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "NPlex";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Specify the plexity of Parent2 used for QTL analysis.";

        /// <summary>
        /// plex数の規定値
        /// </summary>
        public const int DEFAULT = 1;


        /// <summary>
        /// plex数の最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// plex数の仮の最大値
        /// </summary>
        private const int MAXIMUM = 100;


        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// P2 Plex数オプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public Parent2PlexNumberOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _settingValue.Parent2PlexNumber && _settingValue.Parent2PlexNumber <= MAXIMUM)
                return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"Plexity of parent2 should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.Parent2PlexNumber.ToString();

        protected override void SetValue(string value) => _settingValue.Parent2PlexNumber = int.Parse(value);
    }
}
