using PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO;

namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    // 型でサンプルを区別するために継承して型を変えている

    /// <summary>
    /// VCFに記載されているParent1情報
    /// </summary>
    internal class VcfParent1 : VcfSample
    {
        /// <summary>
        /// VCF Parent1情報を作成する。
        /// </summary>
        /// <param name="gt">GT</param>
        /// <param name="ad">AD</param>
        /// <param name="allele">アレル</param>
        public VcfParent1(GT gt, AD ad, string allele)
            : base(gt, ad, allele)
        {
        }
    }
}
