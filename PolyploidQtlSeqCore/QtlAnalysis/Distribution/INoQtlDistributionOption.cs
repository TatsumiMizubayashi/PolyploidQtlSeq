namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// QTLなし分布オプションインターフェイス
    /// </summary>
    public interface INoQtlDistributionOption
    {
        /// <summary>
        /// 倍数性の値を取得する。
        /// </summary>
        int Ploidy { get; }

        /// <summary>
        /// P2のplex数を取得する。
        /// </summary>
        int Parent2PlexNumber { get; }

        /// <summary>
        /// Bulk1の個体数を取得する。
        /// </summary>
        int Bulk1Number { get; }

        /// <summary>
        /// Bulk2の個体数を取得する。
        /// </summary>
        int Bulk2Number { get; }

        /// <summary>
        /// 分布作成時の試行回数を取得する。
        /// </summary>
        int ReplicatesNumber { get; }
    }
}
