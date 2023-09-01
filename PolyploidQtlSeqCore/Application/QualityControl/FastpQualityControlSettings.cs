using PolyploidQtlSeqCore.QualityControl;

namespace PolyploidQtlSeqCore.Application.QualityControl
{
    /// <summary>
    /// Fastp QC設定
    /// </summary>
    internal class FastpQualityControlSettings
    {
        /// <summary>
        ///  Fastp QC設定インスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">Fastp QC設定値</param>
        public FastpQualityControlSettings(IFastpQualityControlSettingValue settingValue)
        {
            InputRawFastqDirectory = new InputRawFastqDirectory(settingValue.InputDir);
            OutputDirectory = new OutputDirectory(settingValue.OutputDir);
            ReadLengthRequired = new ReadLengthRequired(settingValue.ReadLengthRequired);
            NBaseLimit = new NBaseLimit(settingValue.NBaseLimit);
            BaseQuality = new BaseQuality(settingValue.BaseQuality);
            CutTailMeanQuality = new CutTailMeanQuality(settingValue.CutTailMeanQuality);
            CutTailWindowSize = new CutTailWindowSize(settingValue.CutTailWindowSize);
            ThreadNumber = new ThreadNumber(settingValue.ThreadNumber);
        }

        /// <summary>
        /// 入力Fastqディレクトリ
        /// </summary>
        public InputRawFastqDirectory InputRawFastqDirectory { get; }

        /// <summary>
        /// 出力ディレクトリ
        /// </summary>
        public OutputDirectory OutputDirectory { get; }

        /// <summary>
        /// トリム後リード長の最低値
        /// </summary>
        public ReadLengthRequired ReadLengthRequired { get; }

        /// <summary>
        /// N塩基数の上限値
        /// </summary>
        public NBaseLimit NBaseLimit { get; }

        /// <summary>
        /// 塩基クオリティのしきい値
        /// </summary>
        public BaseQuality BaseQuality { get; }

        /// <summary>
        /// 3'末端トリム時の平均クオリティしきい値
        /// </summary>
        public CutTailMeanQuality CutTailMeanQuality { get; }

        /// <summary>
        /// 3'末端トリム時のウインドウサイズ
        /// </summary>
        public CutTailWindowSize CutTailWindowSize { get; }

        /// <summary>
        /// 使用するスレッド数
        /// </summary>
        public ThreadNumber ThreadNumber { get; }


        /// <summary>
        /// FastpCommonOptionに変換する。
        /// </summary>
        /// <returns>FastpCommonOption</returns>
        public FastpCommonOption ToFastpCommonOption()
        {
            return new FastpCommonOption(
                ReadLengthRequired,
                NBaseLimit,
                BaseQuality,
                CutTailMeanQuality,
                CutTailWindowSize,
                OutputDirectory,
                ThreadNumber);
        }
    }
}
