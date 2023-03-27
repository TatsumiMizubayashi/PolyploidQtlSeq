namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    // 型でサンプルを区別するために継承して型を変えている

    /// <summary>
    /// VCFに記載されているParent2情報
    /// </summary>
    internal class VcfParent2 : VcfSample
    {
        /// <summary>
        /// VCF Parent2情報を作成する。
        /// </summary>
        /// <param name="gt">GT値</param>
        /// <param name="allele">アレル</param>
        /// <param name="refCount">Refリード数</param>
        /// <param name="altCount">Altリード数</param>
        public VcfParent2(string gt, string allele, int refCount, int altCount)
            : base(gt, allele, refCount, altCount)
        {
        }
    }
}
