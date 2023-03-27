namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    // 型でサンプルを区別するために継承して型を変えている

    /// <summary>
    /// VCFに記載されているBulk2情報
    /// </summary>
    internal class VcfBulk2 : VcfSample
    {
        /// <summary>
        /// VCF Bulk2情報を作成する。
        /// </summary>
        /// <param name="gt">GT値</param>
        /// <param name="allele">アレル</param>
        /// <param name="refCount">Refリード数</param>
        /// <param name="altCount">Altリード数</param>
        public VcfBulk2(string gt, string allele, int refCount, int altCount)
            : base(gt, allele, refCount, altCount)
        {
        }
    }
}
