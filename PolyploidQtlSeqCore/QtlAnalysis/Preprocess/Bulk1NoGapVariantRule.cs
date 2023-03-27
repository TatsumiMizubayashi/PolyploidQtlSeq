namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    /// <summary>
    /// Bulk1NoGapルール
    /// </summary>
    internal class Bulk1NoGapVariantRule : IAnalyzableVariantRule
    {
        public bool Ok(VcfVariant variant)
        {
            return variant.Bulk1.Depth > 0;
        }
    }
}
