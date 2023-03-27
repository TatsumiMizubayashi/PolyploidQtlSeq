namespace PolyploidQtlSeqCore.QtlAnalysis.Chr
{
    /// <summary>
    /// 解析対象染色体オプションインターフェイス
    /// </summary>
    public interface IAnalysisChrOptions
    {
        /// <summary>
        /// 染色体サイズのしきい値を取得する。
        /// ChrNamesが指定された場合は無視される。
        /// </summary>
        int ChrSizeThreshold { get; }

        /// <summary>
        /// 解析対象染色体名（複数ある場合はカンマ区切り)を取得する。
        /// </summary>
        string AnalysisChrNames { get; }
    }
}
