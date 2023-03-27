namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// Bulk2 F1サンプル
    /// </summary>
    internal class Bulk2F1Samples
    {
        /// <summary>
        /// Bulk2 F1サンプルを作成する。
        /// </summary>
        /// <param name="f1Samples">F1サンプル</param>
        /// <param name="bulk">Bulk2</param>
        public Bulk2F1Samples(F1Samples f1Samples, Bulk2 bulk)
        {
            var reads = f1Samples.ToReads(bulk.Depth);
            SnpIndex = new Bulk2SnpIndex(reads.SnpIndex);
        }

        public Bulk2SnpIndex SnpIndex { get; }
    }
}
