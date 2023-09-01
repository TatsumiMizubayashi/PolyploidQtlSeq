using PolyploidQtlSeqCore.Application.QualityControl;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QualityControl
{
    /// <summary>
    /// 3'末端トリム ウインドウサイズオプション
    /// </summary>
    internal class CutTailWindowSizeOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "W";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "cutTailWindowSize";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Window size when trimmed at 3' end. cut_tail_window_size in fastp.";


        /// <summary>
        /// 3'末端トリム時のウインドウサイズの規定値
        /// </summary>
        public const int DEFAULT = 1;

        /// <summary>
        /// 3'末端トリム時のウインドウサイズの最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 3'末端トリム時のウインドウサイズの最大値
        /// </summary>
        private const int MAXIMUM = 100;

        private readonly IFastpQualityControlSettingValue _optionValue;

        /// <summary>
        /// 3'末端トリム ウインドウサイズオプションインスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">Fastq QCオプション値</param>
        public CutTailWindowSizeOption(IFastpQualityControlSettingValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _optionValue.CutTailWindowSize && _optionValue.CutTailWindowSize <= MAXIMUM) return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME,
                $"Window size when trimmed at 3' end should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.CutTailWindowSize.ToString();

        protected override void SetValue(string value) => _optionValue.CutTailWindowSize = int.Parse(value);
    }
}
