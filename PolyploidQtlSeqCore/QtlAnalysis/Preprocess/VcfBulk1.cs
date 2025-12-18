using PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO;

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
        /// <param name="gt">GT</param>
        /// <param name="ad">AD</param>
        /// <param name="allele">アレル</param>
        public VcfBulk1(GT gt, AD ad, string allele)
            : base(gt, ad, allele)
        {
        }
    }
}
