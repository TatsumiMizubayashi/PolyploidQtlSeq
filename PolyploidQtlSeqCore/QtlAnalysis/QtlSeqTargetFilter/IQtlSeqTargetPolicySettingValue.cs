namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// QTL seq対象ポリシー設定値のインターフェイス
    /// </summary>
    public interface IQtlSeqTargetPolicySettingValue
    {
        /// <summary>
        /// P1 MostAlleleRateのしきい値を取得する。
        /// </summary>
        double Parent1MostAlleleRateThreshold { get; }

        /// <summary>
        /// P2 SNP-index範囲を取得する。
        /// </summary>
        string Parent2SnpIndexRange { get; }

        /// <summary>
        /// 最低Depthのしきい値を取得する。
        /// </summary>
        int MinimumDepthThreshold { get; }

        /// <summary>
        /// 最大Bulk SNP-indexのしきい値を取得する。
        /// </summary>
        double MaxBulkSnpIndexThreshold { get; }
    }
}
