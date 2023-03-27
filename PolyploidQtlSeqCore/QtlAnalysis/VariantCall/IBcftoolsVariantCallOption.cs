namespace PolyploidQtlSeqCore.QtlAnalysis.VariantCall
{
    /// <summary>
    /// Bcftools変異検出オプションインターフェイス
    /// </summary>
    public interface IBcftoolsVariantCallOption
    {
        /// <summary>
        /// MQ最低値を取得する。
        /// </summary>
        int MinMq { get; }

        /// <summary>
        /// BQ最低値を取得する。
        /// </summary>
        int MinBq { get; }

        /// <summary>
        /// Adjust MQを取得する。
        /// </summary>
        int AdjustMq { get; }
    }
}
