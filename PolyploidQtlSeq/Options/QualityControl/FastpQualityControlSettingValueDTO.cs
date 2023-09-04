using PolyploidQtlSeqCore.Application.QualityControl;

namespace PolyploidQtlSeq.Options.QualityControl
{
    /// <summary>
    /// Fastp QC設定値DTO
    /// </summary>
    internal class FastpQualityControlSettingValueDTO : IFastpQualityControlSettingValue
    {
        /// <summary>
        /// Fastp QC設定値DTOインスタンスを作成する。
        /// </summary>
        public FastpQualityControlSettingValueDTO()
        {
            InputDir = "";
            OutputDir = "";
            ReadLengthRequired = 0;
            NBaseLimit = 0;
            BaseQuality = 0;
            CutTailMeanQuality = 0;
            CutTailWindowSize = 0;
            ThreadNumber = 0;
        }

        public string InputDir { get; init; }

        public string OutputDir { get; init; }

        public int ReadLengthRequired { get; init; }

        public int NBaseLimit { get; init; }

        public int BaseQuality { get; init; }

        public int CutTailMeanQuality { get; init; }

        public int CutTailWindowSize { get; init; }

        public int ThreadNumber { get; init; }
    }
}
