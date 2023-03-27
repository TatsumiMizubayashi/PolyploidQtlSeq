namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess
{
    /// <summary>
    /// 親1GTホモ型ルール
    /// </summary>
    internal class Parent1GtHomoRule : IAnalyzableVariantRule
    {
        public bool Ok(VcfVariant variant)
        {
            return variant.Parent1.GT == GtType.RefHomo
                || variant.Parent1.GT == GtType.AltHomo;
        }
    }
}
