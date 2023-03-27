namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    /// <summary>
    /// 親1NoGapルール
    /// </summary>
    internal class Parent1NoGapVariantRule : IAnalyzableVariantRule
    {
        public bool Ok(VcfVariant variant)
        {
            return variant.Parent1.Depth > 0;
        }
    }
}
