using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// 分布作成試行回数オプション
    /// </summary>
    internal class ReplicatesNumberOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "N";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "NRep";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Number of simulation replicates to generate a null distribution which is free from QTLs.";

        /// <summary>
        /// 試行回数の規定値
        /// </summary>
        public const int DEFAULT = 5000;

        /// <summary>
        /// 試行回数の最小値
        /// </summary>
        private const int MINIMUM = 1000;

        /// <summary>
        /// 試行回数の最大値
        /// </summary>
        private const int MAXIMUM = 1_000_000;

        private readonly IQtlSeqAnalysisOptionValue _settingValue;

        /// <summary>
        /// 分布作成試行回数オプションインスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">QTL-Seq解析設定値</param>
        public ReplicatesNumberOption(IQtlSeqAnalysisOptionValue settingValue)
        {
            _settingValue = settingValue;
        }   

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _settingValue.ReplicatesNumber && _settingValue.ReplicatesNumber <= MAXIMUM)
                return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"Number of simulation replicates should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _settingValue.ReplicatesNumber.ToString();

        protected override void SetValue(string value) => _settingValue.ReplicatesNumber = int.Parse(value);
    }
}
