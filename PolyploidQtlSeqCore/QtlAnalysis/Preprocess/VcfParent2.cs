using PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO;

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
        /// <param name="gt">GT</param>
        /// <param name="ad">AD</param>
        /// <param name="allele">アレル</param>
        public VcfParent2(GT gt, AD ad, string allele)
            : base(gt, ad, allele)
        {
        }
    }
}
