using PolyploidQtlSeqCore.QtlAnalysis.Preprocess;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// SNP-index計算機
    /// </summary>
    internal static class SnpIndexCalculator
    {
        /*
         * SNP-index = P2型リード数 / Depth
         */

        /// <summary>
        /// サンプルのSNP-indexを計算する。
        /// </summary>
        /// <param name="vcfP1">Parent1</param>
        /// <param name="vcfSample">サンプル</param>
        /// <returns>SNP-index</returns>
        public static double Calc(VcfParent1 vcfP1, VcfSample vcfSample)
        {
            return vcfP1.GT == GtType.RefHomo
                ? CalcParent1RefHomoSnpIndex(vcfSample)
                : vcfP1.GT == GtType.AltHomo
                    ? CalcParent1AltHomoSnpIndex(vcfSample)
                    : throw new InvalidOperationException("P1 is not homotypic.");
        }

        private static double CalcParent1RefHomoSnpIndex(VcfSample sample)
        {
            // P1がRefホモ型なのでP2はAlt型
            return sample.AltCount / (double)sample.Depth;
        }

        private static double CalcParent1AltHomoSnpIndex(VcfSample sample)
        {
            // P1がAltホモ型なのでP2はRef型
            return sample.RefCount / (double)sample.Depth;
        }
    }
}
