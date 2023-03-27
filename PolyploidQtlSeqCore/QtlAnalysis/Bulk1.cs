using PolyploidQtlSeqCore.QtlAnalysis.Preprocess;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// Bulk1
    /// </summary>
    internal class Bulk1
    {
        /// <summary>
        /// Bulk1を作成する。
        /// </summary>
        /// <param name="vcfP1">VcfP1</param>
        /// <param name="vcfBulk1">VcfBulk1</param>
        public Bulk1(VcfParent1 vcfP1, VcfBulk1 vcfBulk1)
        {
            Depth = vcfBulk1.Depth;
            Allele = vcfBulk1.Allele;

            var snpIndexValue = SnpIndexCalculator.Calc(vcfP1, vcfBulk1);
            SnpIndex = new Bulk1SnpIndex(snpIndexValue);
        }

        /// <summary>
        /// SNP-indexを取得する。
        /// </summary>
        public Bulk1SnpIndex SnpIndex { get; }

        /// <summary>
        /// Depthを取得する。
        /// </summary>
        public int Depth { get; }

        /// <summary>
        /// アレル塩基を取得する。
        /// </summary>
        public string Allele { get; }
    }
}
