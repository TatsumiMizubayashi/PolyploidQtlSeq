namespace PolyploidQtlSeqCore.QtlAnalysis.Chr
{
    /// <summary>
    /// 解析対象染色体設定
    /// </summary>
    internal class AnalysisChrSettings
    {
        /// <summary>
        /// 解析対象染色体設定を作成する。
        /// </summary>
        /// <param name="settingValue">設定値</param>
        public AnalysisChrSettings(IAnalysisChrSettingValue settingValue)
        {
            ChrSizeThreshold = new ChrSizeThreshold(settingValue.ChrSizeThreshold);
            AnalysisChrNames = new AnalysisChrNames(settingValue.AnalysisChrNames);
        }

        /// <summary>
        /// 染色体サイズしきい値
        /// </summary>
        public ChrSizeThreshold ChrSizeThreshold { get; }

        /// <summary>
        /// 解析対象染色体名
        /// </summary>
        public AnalysisChrNames AnalysisChrNames { get; }
    }
}
