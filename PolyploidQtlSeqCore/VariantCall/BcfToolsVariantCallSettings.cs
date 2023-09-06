using PolyploidQtlSeqCore.QtlAnalysis;
using PolyploidQtlSeqCore.Share;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// Bcftools 変異検出オプション
    /// </summary>
    internal class BcfToolsVariantCallSettings
    {
        /// <summary>
        /// Bcftools変異検出オプションを作成する。
        /// </summary>
        /// <param name="settingValue">bcftools変異検出設定値</param>
        public BcfToolsVariantCallSettings(IBcftoolsVariantCallSettingValue settingValue)
        {
            ReferenceSequence = new ReferenceSequence(settingValue.ReferenceSequence);
            MinmumBaseQuality = new MinmumBaseQuality(settingValue.MinBq);
            MinmumMappingQuality = new MinmumMappingQuality(settingValue.MinMq);
            AdjustMappingQuality = new AdjustMappingQuality(settingValue.AdjustMq);
            OutputDirectory = new OutputDirectory(settingValue.OutputDir);
            ThreadNumber = new ThreadNumber(settingValue.ThreadNumber);
        }

        /// <summary>
        /// リファレンスシークエンス
        /// </summary>
        public ReferenceSequence ReferenceSequence { get; }

        /// <summary>
        /// Base Quality最低値
        /// </summary>
        public MinmumBaseQuality MinmumBaseQuality { get; }

        /// <summary>
        /// Mapping Quality最低値
        /// </summary>
        public MinmumMappingQuality MinmumMappingQuality { get; }

        /// <summary>
        /// Adjust Mapping Quality
        /// </summary>
        public AdjustMappingQuality AdjustMappingQuality { get; }

        /// <summary>
        /// 出力ディレクトリ
        /// </summary>
        public OutputDirectory OutputDirectory { get; }

        /// <summary>
        /// スレッド数
        /// </summary>
        public ThreadNumber ThreadNumber { get; }
    }
}
