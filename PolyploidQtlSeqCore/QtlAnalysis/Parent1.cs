using PolyploidQtlSeqCore.QtlAnalysis.Preprocess;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// 親1
    /// </summary>
    internal class Parent1
    {
        /// <summary>
        /// Parent1を作成する。
        /// </summary>
        /// <param name="vcfP1">vcfP1</param>
        public Parent1(VcfParent1 vcfP1)
        {
            Allele = vcfP1.Allele;
            Depth = vcfP1.Depth;
            MaxAlleleRate = Math.Max(vcfP1.RefCount, vcfP1.AltCount) / (double)vcfP1.Depth;
        }

        /// <summary>
        /// アレル塩基を取得する。
        /// </summary>
        public string Allele { get; }

        /// <summary>
        /// Depthを取得する。
        /// </summary>
        public int Depth { get; }

        /// <summary>
        /// 最多アレル割合を取得する。
        /// </summary>
        public double MaxAlleleRate { get; }
    }
}
