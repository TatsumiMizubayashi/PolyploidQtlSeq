namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// 親1 MostAllele割合ルール
    /// </summary>
    internal class Parent1MostAlleleRateRule : IQtlSeqTargetVariantRule
    {
        private readonly Parent1MostAlleleRateThreshold _threshold;

        /// <summary>
        /// 親1 MostAllele割合ルールを作成する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        public Parent1MostAlleleRateRule(Parent1MostAlleleRateThreshold threshold)
        {
            _threshold = threshold;
        }


        public bool Ok(SnpIndexVariant variant)
        {
            return variant.Parent1.MaxAlleleRate >= _threshold.Value;
        }
    }
}
