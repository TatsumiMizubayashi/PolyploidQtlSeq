namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// QTL判定しきい値のΔSNP-index
    /// </summary>
    internal class QtlThresholdDeltaSnpIndex
    {
        /// <summary>
        /// QTL判定しきい値のΔSNP-indexを作成する。
        /// </summary>
        /// <param name="threshold">ΔSNP-indexしきい値</param>
        public QtlThresholdDeltaSnpIndex(double threshold)
        {
            if (threshold < 0) throw new ArgumentException(null, nameof(threshold));

            Value = threshold;
        }

        /// <summary>
        /// しきい値のΔSNP-indexを取得する。
        /// </summary>
        public double Value { get; }
    }
}
