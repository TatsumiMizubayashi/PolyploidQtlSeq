namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    /// <summary>
    /// Bulk1NoGapルール
    /// </summary>
    internal class Bulk2NoGapVariantRule : IAnalyzableVariantRule
    {
        public bool Ok(VcfVariant variant)
        {
            return variant.Bulk2.Depth > 0;
        }
    }
}
