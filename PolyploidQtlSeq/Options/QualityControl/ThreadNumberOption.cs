using PolyploidQtlSeqCore.Application.QualityControl;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QualityControl
{
    /// <summary>
    /// スレッド数オプション
    /// </summary>
    internal class ThreadNumberOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "t";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "thread";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Number of threads to use. Up to 16.";

        /// <summary>
        /// 使用するスレッド数の規定値
        /// </summary>
        public const int DEFAULT = 10;

        /// <summary>
        /// 使用するスレッド数の最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 使用するスレッド数の最大値
        /// </summary>
        private const int MAXIMUM = 16;

        private readonly IFastpQualityControlOptionValue _optionValue;

        /// <summary>
        /// スレッド数オプション インスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">Fastq QCオプション値</param>
        public ThreadNumberOption(IFastpQualityControlOptionValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _optionValue.ThreadNumber && _optionValue.ThreadNumber <= MAXIMUM) return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"Number of threads should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.ThreadNumber.ToString();

        protected override void SetValue(string value) => _optionValue.ThreadNumber = int.Parse(value);
    }
}
