using PolyploidQtlSeqCore.QtlAnalysis.Preprocess;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// 親2
    /// </summary>
    internal class Parent2
    {
        /// <summary>
        /// Parent2を作成する。
        /// </summary>
        /// <param name="vcfP1">VcfP1</param>
        /// <param name="vcfP2">VcfP2</param>
        public Parent2(VcfParent1 vcfP1, VcfParent2 vcfP2)
        {
            Depth = vcfP2.Depth;
            SnpIndex = SnpIndexCalculator.Calc(vcfP1, vcfP2);
        }

        /// <summary>
        /// SNP-indexを取得する。
        /// </summary>
        public double SnpIndex { get; }

        /// <summary>
        /// Depthを取得する。
        /// </summary>
        public int Depth { get; }
    }
}
