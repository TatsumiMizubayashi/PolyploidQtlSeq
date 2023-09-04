using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// 倍数性オプション
    /// </summary>
    internal class PloidyOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "p";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "ploidy";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Ploidy.";

        /// <summary>
        /// ploidyの規定値
        /// </summary>
        public const int DEFAULT = 4;

        /// <summary>
        /// ploidyの最小値
        /// </summary>
        private const int MINIMUM = 2;

        /// <summary>
        /// ploidyの最大値
        /// </summary>
        private const int MAXIMUM = 20;


        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// 倍数性オプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public PloidyOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }

        public override DataValidationResult Validation()
        {
            if(_settingValue.Ploidy < MINIMUM || _settingValue.Ploidy > MAXIMUM)
                return new DataValidationResult(SHORT_NAME, LONG_NAME, $"Ploidy should be an integer between {MINIMUM} and {MAXIMUM}.");

            if (_settingValue.Ploidy % 2 != 0)
                return new DataValidationResult(SHORT_NAME, LONG_NAME, "Specify an even number for Ploidy.");

            return new DataValidationResult();
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.Ploidy.ToString();

        protected override void SetValue(string value) => _settingValue.Ploidy = int.Parse(value);
    }
}
