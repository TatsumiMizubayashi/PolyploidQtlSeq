namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// Bulk2 SNP-index
    /// </summary>
    internal class Bulk2SnpIndex
    {
        /// <summary>
        /// Bulk2 SNP-indexを作成する。
        /// </summary>
        /// <param name="snpIndex">SNP-index</param>
        public Bulk2SnpIndex(double snpIndex)
        {
            Value = snpIndex;
        }

        /// <summary>
        /// SNP-indexを取得する。
        /// </summary>
        public double Value { get; }

        /// <summary>
        /// ΔSNP-indexを計算する。
        /// </summary>
        /// <param name="bulk1SnpIndex">Bulk1 SNP-index</param>
        /// <returns>ΔSNP-index</returns>
        public DeltaSnpIndex CalcDeltaSnpIndex(Bulk1SnpIndex bulk1SnpIndex)
        {
            return new DeltaSnpIndex(Value - bulk1SnpIndex.Value);
        }
    }
}
