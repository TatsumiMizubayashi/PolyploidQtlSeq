namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// Bulk1 F1サンプル
    /// </summary>
    internal class Bulk1F1Samples
    {
        /// <summary>
        /// Bulk1 F1サンプルを作成する。
        /// </summary>
        /// <param name="f1Samples">F1サンプル</param>
        /// <param name="bulk">Bulk1情報</param>
        public Bulk1F1Samples(F1Samples f1Samples, Bulk1 bulk)
        {
            var reads = f1Samples.ToReads(bulk.Depth);
            SnpIndex = new Bulk1SnpIndex(reads.SnpIndex);
        }

        /// <summary>
        /// Bulk1 SNP-indexを取得する。
        /// </summary>
        public Bulk1SnpIndex SnpIndex { get; }
    }
}
