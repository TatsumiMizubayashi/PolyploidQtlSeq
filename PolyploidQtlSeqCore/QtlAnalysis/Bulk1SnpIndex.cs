namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// Bulk1SNP-index
    /// </summary>
    internal class Bulk1SnpIndex
    {
        /// <summary>
        /// Bulk1SnpIndexを作成する。
        /// </summary>
        /// <param name="snpIndex">SNP-index</param>
        public Bulk1SnpIndex(double snpIndex)
        {
            Value = snpIndex;
        }

        /// <summary>
        /// SNP-indexを取得する。
        /// </summary>
        public double Value { get; }
    }
}
