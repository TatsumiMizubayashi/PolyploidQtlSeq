using PolyploidQtlSeqCore.Application.QualityControl;
using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeq.Options.QualityControl
{
    /// <summary>
    /// 必須リード長オプション
    /// </summary>
    internal class ReadLengthRequiredOption : Option
    {
        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "l";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "lengthRequired";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Minimum read length after trimming. length_required in fastp.";

        /// <summary>
        /// リード最低長の規定値
        /// </summary>
        public const int DEFAULT = 50;

        /// <summary>
        /// リード最低長の最小値
        /// </summary>
        private const int MINIMUM = 10;

        /// <summary>
        /// リード最低長の最大値
        /// </summary>
        private const int MAXIMUM = 300;

        private readonly IFastpQualityControlSettingValue _optionValue;

        /// <summary>
        /// 必須リード長オプション インスタンスを作成する。
        /// </summary>
        /// <param name="optionValue">Fastq QCオプション値</param>
        public ReadLengthRequiredOption(IFastpQualityControlSettingValue optionValue)
        {
            _optionValue = optionValue;
        }

        public override DataValidationResult Validation()
        {
            if (MINIMUM <= _optionValue.ReadLengthRequired && _optionValue.ReadLengthRequired <= MAXIMUM) return new DataValidationResult();

            return new DataValidationResult(SHORT_NAME, LONG_NAME, 
                $"Minimum read length should be an integer between {MINIMUM} and {MAXIMUM}.");
        }

        protected override string GetLongName() => LONG_NAME;

        protected override string GetShortName() => SHORT_NAME;

        protected override string GetStringValue() => _optionValue.ReadLengthRequired.ToString();

        protected override void SetValue(string value) => _optionValue.ReadLengthRequired = int.Parse(value);
    }
}
