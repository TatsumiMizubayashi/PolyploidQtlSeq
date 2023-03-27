namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    // 型でサンプルを区別するために継承して型を変えている

    /// <summary>
    /// VCFに記載されているBulk1情報
    /// </summary>
    internal class VcfBulk1 : VcfSample
    {
        /// <summary>
        /// VCF Bulk1情報を作成する。
        /// </summary>
        /// <param name="gt">GT値</param>
        /// <param name="allele">アレル</param>
        /// <param name="refCount">Refリード数</param>
        /// <param name="altCount">Altリード数</param>
        public VcfBulk1(string gt, string allele, int refCount, int altCount)
            : base(gt, allele, refCount, altCount)
        {
        }
    }
}
