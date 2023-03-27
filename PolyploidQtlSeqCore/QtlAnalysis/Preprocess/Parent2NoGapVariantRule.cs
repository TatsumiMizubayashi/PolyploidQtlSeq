namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    /// <summary>
    /// 親2NoGapルール
    /// </summary>
    internal class Parent2NoGapVariantRule : IAnalyzableVariantRule
    {
        public bool Ok(VcfVariant variant)
        {
            return variant.Parent2.Depth > 0;
        }
    }
}
