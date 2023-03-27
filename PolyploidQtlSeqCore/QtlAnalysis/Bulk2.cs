using PolyploidQtlSeqCore.QtlAnalysis.Preprocess;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// Bulk2
    /// </summary>
    internal class Bulk2
    {
        /// <summary>
        /// Bulk2を作成する。
        /// </summary>
        /// <param name="vcfP1">VcfP1</param>
        /// <param name="vcfBulk2">VcfP2</param>
        public Bulk2(VcfParent1 vcfP1, VcfBulk2 vcfBulk2)
        {
            Depth = vcfBulk2.Depth;
            Allele = vcfBulk2.Allele;

            var snpIndexValue = SnpIndexCalculator.Calc(vcfP1, vcfBulk2);
            SnpIndex = new Bulk2SnpIndex(snpIndexValue);
        }

        /// <summary>
        /// SNP-indexを取得する。
        /// </summary>
        public Bulk2SnpIndex SnpIndex { get; }

        /// <summary>
        /// Depthを取得する。
        /// </summary>
        public int Depth { get; }

        /// <summary>
        /// アレル塩基を取得する。
        /// </summary>
        public string Allele { get; }

        /// <summary>
        /// ΔSNP-indexを計算する。
        /// </summary>
        /// <param name="bulk1">Bulk1</param>
        /// <returns>ΔSNP-index</returns>
        public DeltaSnpIndex CalcDeltaSnpIndex(Bulk1 bulk1)
        {
            return SnpIndex.CalcDeltaSnpIndex(bulk1.SnpIndex);
        }
    }
}
