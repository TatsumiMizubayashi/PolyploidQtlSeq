using PolyploidQtlSeqCore.QtlAnalysis.Chr;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// 変異検出シナリオ設定
    /// </summary>
    internal class VariantCallScenarioSettings
    {
        /// <summary>
        /// 変異検出シナリオ設定インスタンスを作成する。
        /// </summary>
        /// <param name="settingValue">変異検出シナリオ設定値</param>
        public VariantCallScenarioSettings(IVariantCallScenarioSettingValue settingValue)
        {
            AnalysisChrSettings = new AnalysisChrSettings(settingValue);
            BcfToolsVariantCallSettings = new BcfToolsVariantCallSettings(settingValue);
            SnpEffSettings = new SnpEffSettings(settingValue);
        }

        /// <summary>
        /// 解析染色体設定を取得する。
        /// </summary>
        public AnalysisChrSettings AnalysisChrSettings { get; }

        /// <summary>
        /// bcftools変異検出設定を取得する。
        /// </summary>
        public BcfToolsVariantCallSettings BcfToolsVariantCallSettings { get; }

        /// <summary>
        /// SnpEff設定を取得する。
        /// </summary>
        public SnpEffSettings SnpEffSettings { get; }
    }
}
